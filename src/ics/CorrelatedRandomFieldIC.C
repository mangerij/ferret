/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/



// Credits to A. Hagerstrom (NIST) for this piece of C++ code 

//NOTE: this only works for a regular geometry and not within subdimensional boundary restricted computational domains.


#include "CorrelatedRandomFieldIC.h"


#include "SubProblem.h" //need both of these to pull total node number
#include "MooseMesh.h" 

#include "libmesh/point.h"
#include <cmath>

#include "MooseRandom.h"

#include "libmesh/point.h"

/**
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
**/

namespace
{
/**
 * Method for retrieving or generating and caching a value in a map.
 */
inline Real
valueHelper(dof_id_type id, MooseRandom & generator, std::map<dof_id_type, Real> & map)
{
  auto it_pair = map.lower_bound(id);

  // Do we need to generate a new number?
  if (it_pair == map.end() || it_pair->first != id)
    it_pair = map.emplace_hint(it_pair, id, generator.rand(id));

  return it_pair->second;
}
}


registerMooseObject("FerretApp", CorrelatedRandomFieldIC);

template<>
InputParameters validParams<CorrelatedRandomFieldIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("Lcorr", "correlation length scale");

  MooseEnum dims("1=1 2 3");
  params.addRequiredParam<MooseEnum>("dim", dims, "The dimension of the mesh to recieve correlated field"); 
  params.addParam<Real>("xmin", 0.0, "Lower X Coordinate of the generated mesh");
  params.addParam<Real>("ymin", 0.0, "Lower Y Coordinate of the generated mesh");
  params.addParam<Real>("zmin", 0.0, "Lower Z Coordinate of the generated mesh");
  params.addParam<Real>("xmax", 1.0, "Upper X Coordinate of the generated mesh");
  params.addParam<Real>("ymax", 1.0, "Upper Y Coordinate of the generated mesh");
  params.addParam<Real>("zmax", 1.0, "Upper Z Coordinate of the generated mesh");
  params.addParam<Real>("Nnodes", 1.0, "total number of nodes to loop over");
  params.addParam<unsigned int>("seed", 0, "Seed value for the random number generator"); 
  return params;
}

CorrelatedRandomFieldIC::CorrelatedRandomFieldIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    _is_nodal(_var.isNodal()),
    _use_legacy(getParam<bool>("legacy_generator")),
    _elem_random_generator(nullptr),
    _node_random_generator(nullptr),
    _Lcorr(getParam<Real>("Lcorr")),
    _dim(getParam<MooseEnum>("dim")),
    _xmin(getParam<Real>("xmin")),
    _xmax(getParam<Real>("xmax")),
    _ymin(getParam<Real>("ymin")),
    _ymax(getParam<Real>("ymax")),
    _zmin(getParam<Real>("zmin")),
    _zmax(getParam<Real>("zmax")),
    _Nnodes(getParam<Real>("Nnodes")),
    _mesh(_fe_problem.mesh().getMesh())
{
  unsigned int processor_seed = getParam<unsigned int>("seed");
  MooseRandom::seed(processor_seed);

  if (_use_legacy)
  {
    auto proc_id = processor_id();
    if (proc_id > 0)
    {
      for (processor_id_type i = 0; i < proc_id; ++i)
        processor_seed = MooseRandom::randl();
      MooseRandom::seed(processor_seed);
    }
  }
  else
  {
    _elem_random_data =
        libmesh_make_unique<RandomData>(_fe_problem, false, EXEC_INITIAL, MooseRandom::randl());
    _node_random_data =
        libmesh_make_unique<RandomData>(_fe_problem, true, EXEC_INITIAL, MooseRandom::randl());

    _elem_random_generator = &_elem_random_data->getGenerator();
    _node_random_generator = &_node_random_data->getGenerator();
  }
}


void
CorrelatedRandomFieldIC::initialSetup()
{
  if (!_use_legacy)
  {
    _elem_random_data->updateSeeds(EXEC_INITIAL);
    _node_random_data->updateSeeds(EXEC_INITIAL);
  }
}

Real
CorrelatedRandomFieldIC::generateRandom()
{
  Real rand_num;

  if (_use_legacy)
  {
    mooseDeprecated("legacy_generator is deprecated. Please set \"legacy_generator = false\". This "
                    "capability will be removed after 11/01/2018");

    // Random number between 0 and 1
    rand_num = MooseRandom::rand();
  }
  else
  {
    if (_current_node)
      rand_num = valueHelper(_current_node->id(), *_node_random_generator, _node_numbers);
    else if (_current_elem)
      rand_num = valueHelper(_current_elem->id(), *_elem_random_generator, _elem_numbers);
    else
      mooseError("We can't generate parallel consistent random numbers for this kind of variable "
                 "yet. Please contact the MOOSE team for assistance");
  }

  return rand_num;
}


//stores the coefficient matrix
std::vector<std::vector<std::vector<std::vector<Real>>>>
CorrelatedRandomFieldIC::fourierCoeffs()
{
    int nodes = _Nnodes; //_mesh.n_nodes();
    std::cout << nodes;
    std::array<Real, LIBMESH_DIM> L = {{_xmax - _xmin, _dim > 1 ? _ymax - _ymin : 0, _dim > 2 ? _zmax - _zmin : 0}};
  
   std::vector<std::vector<std::vector<std::vector<Real>>>> output = std::vector<std::vector<std::vector<std::vector<Real>>>>(2*nodes+1,std::vector<std::vector<std::vector<Real>>>(2*nodes+1,std::vector<std::vector<Real>>(2*nodes+1,std::vector<Real>(2,0))));
    
    Real rand1  = 0.0;
    Real rand2  = 0.0;
    
    Real nrand1  = 0.0;
    Real nrand2  = 0.0;   //random seed needs to be MOOSE-y but for now we can test this
    
    Real randnorm = 0.0;
    Real amp = 0.0;
    
    RealVectorValue kv(0.0, 0.0, 0.0); //kx ky kz
    Real k2 = 0.0;

    Real norm = std::sqrt(8.0*libMesh::pi/(_Lcorr*L[0]*L[1]*L[2]));


   // srand(seed);
    for (unsigned int i = 0; i<2*nodes+1; i++){
        for (unsigned int j = 0;j<2*nodes+1;j++){
            for (unsigned int k = 0; k<2*nodes+1; k++){
                // the fourier transform of a real variable 
                // has the symmetry F(i,j,k) = F^*(-i,-j,-k)
                
                if (i+j+k < 0){
                    continue;
                }
                
                kv(0) = 2.0*libMesh::pi*(i-nodes)/L[0];
                kv(1) = 2.0*libMesh::pi*(j-nodes)/L[1];
                kv(2) = 2.0*libMesh::pi*(k-nodes)/L[2];
                
                // compute absolute value of coefficient
                k2 = kv(0)*kv(0) + kv(1)*kv(1) + kv(2)*kv(2);
                amp = norm/(k2 + std::pow(1.0/_Lcorr,2));
                
                // Box-Muller method for gaussian random numbers
                
                rand1  = (generateRandom()+1)/1000000.0;
                //Moose::out << "\n ";
                //Moose::out << "\n rand1 = "; std::cout << rand1;
                //Moose::out << "\n ";
                rand2  = (generateRandom()+1)/1000000.0;
                
                //construct random numbers with unit variance
                //with unit variance for nrand1^2 + nrand2^2
                nrand1 = sqrt(-2.0*std::log(rand1))*cos(2.0*libMesh::pi*rand2)/sqrt(2.0);

                //Moose::out << "\n nrand1 = "; std::cout << nrand1;

                nrand2 = sqrt(-2.0*std::log(rand1))*sin(2.0*libMesh::pi*rand2)/sqrt(2.0);

                //Moose::out << "\n nrand2 = "; std::cout << nrand2;

                //enforce 0 mean
                if (i == nodes && j == nodes && k == nodes){
                    amp = 0.0;
                }

                output[i][j][k][0] = amp*nrand1;
                output[i][j][k][1] = amp*nrand2;
            }
        }
    }
  //  Moose::out << "\n ";
   // Moose::out << "\n output[0][0][0][0] = "; std::cout << output[0][0][0][0];
  //  Moose::out << "\n ";
   // Moose::out << "\n ";
  //  Moose::out << "\n output[0][0][0][1] = "; std::cout << output[0][0][0][1];
   // Moose::out << "\n ";
    return output;
}

Real
CorrelatedRandomFieldIC::evaluate(const Point & p, std::vector<std::vector<std::vector<std::vector<Real>>>> output)// array4ptr coeffs???) //this is main()
{
  int nodes = _mesh.n_nodes();
  Real result = 0.0;
  RealVectorValue kv(0.0, 0.0, 0.0); //kx ky kz
  Real phase = 0.0;

  std::array<Real, LIBMESH_DIM> L = {{_xmax - _xmin, _dim > 1 ? _ymax - _ymin : 0, _dim > 2 ? _zmax - _zmin : 0}};

  for (unsigned int i = 0; i<2*nodes+1; i++)
  {
    for (unsigned int j = 0;j<2*nodes+1;j++)
    {
      for (unsigned int k = 0; k<2*nodes+1; k++)
      {
        kv(0) = 2.0*libMesh::pi*double(i-nodes)/L[0];
        kv(1) = 2.0*libMesh::pi*double(j-nodes)/L[1];
        kv(2) = 2.0*libMesh::pi*double(k-nodes)/L[2];
                
        phase = kv(0)*p(0) + kv(1)*p(1) + kv(2)*p(2);
        //the real part
        result += output[i][j][k][0]*cos(phase);
        result += output[i][j][k][1]*sin(phase); //???
        }
     }
  }
  return result;
}



Real
CorrelatedRandomFieldIC::evaluateNoArray(const Point & p)
{
    int nodes = _mesh.n_nodes();
    Real result = 0.0;
    std::array<Real, LIBMESH_DIM> L = {{_xmax - _xmin, _dim > 1 ? _ymax - _ymin : 0, _dim > 2 ? _zmax - _zmin : 0}};

    RealVectorValue kv(0.0, 0.0, 0.0); //kx ky kz

    Real k2 = 0.0;
    Real phase = 0.0;
      
    Real rand1  = 0.0;
    Real rand2  = 0.0;
    
    Real nrand1  = 0.0;
    Real nrand2  = 0.0;
    
    Real randnorm = 0.0;
    Real amp = 0.0;
    
    // normalization constant
    double norm = sqrt(8.0*libMesh::pi/(_Lcorr*L[0]*L[1]*L[2]));
    
    // sum over kspace
    for (int i = 0; i<2*nodes+1; i++)
    {
        for (int j = 0;j<2*nodes+1;j++)
        {
            for (int k = 0; k<2*nodes+1; k++)
            {
                // k vector 
                kv(0) = 2.0*libMesh::pi*double(i-nodes)/L[0];
                kv(1) = 2.0*libMesh::pi*double(j-nodes)/L[1];
                kv(2) = 2.0*libMesh::pi*double(k-nodes)/L[2];
                
                // compute absolute value of coefficient
                k2 = kv(0)*kv(0) + kv(1)*kv(1) + kv(2)*kv(2);
                amp = norm/(k2 + std::pow(1.0/_Lcorr,2));
                
                // Box-Muller method for gaussian random numbers
                rand1  = (generateRandom()+1)/1000000.0;
                rand2  = (generateRandom()+1)/1000000.0;
                
                //construct random numbers with unit variance
                //with unit variance for nrand1^2 + nrand2^2
                nrand1 = sqrt(-2.0*std::log(rand1))*cos(2.0*libMesh::pi*rand2)/sqrt(2.0);
                nrand2 = sqrt(-2.0*std::log(rand1))*sin(2.0*libMesh::pi*rand2)/sqrt(2.0);

                
                //enforce 0 mean
                if (i == nodes && j == nodes && k == nodes){
                    amp = 0.0;
                }
                phase = kv(0)*p(0) + kv(1)*p(1) + kv(2)*p(2);
                
                //the real part
                result+=amp*(nrand1*cos(phase)+nrand2*sin(phase));
            }
        }
    }
    return result;
}

Real
CorrelatedRandomFieldIC::value(const Point & p) //this is main()
{
  //CorrelatedRandomFieldIC::fourierCoeffs();
//  std::cout << d[0][0][0][0];  //issues with this. Why is it -nan??
  return 0.0; //evaluate(p,d);
}
