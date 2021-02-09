/*
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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

   // Credits to A. Hagerstrom (NIST) for this piece of C++ code

**/

#include "CorrelatedRandomFieldIC.h"

#include "SubProblem.h" //need both of these to pull total node number
#include "MooseMesh.h"

#include "libmesh/point.h"
#include <cmath>

#include "MooseRandom.h"

#include "libmesh/point.h"
#include "libmesh/utility.h"
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
  params.addParam<Real>("xmax", 0.0, "Upper X Coordinate of the generated mesh");
  params.addParam<Real>("ymax", 0.0, "Upper Y Coordinate of the generated mesh");
  params.addParam<Real>("zmax", 0.0, "Upper Z Coordinate of the generated mesh");
  params.addParam<Real>("Nnodes", 1.0, "total number of nodes to loop over");
  params.addParam<unsigned int>("seed", 0, "Seed value for the random number generator");
  params.addParam<bool>(
      "legacy_generator",
      false,
      "Determines whether or not the legacy generator (deprecated) should be used.");
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
  dcoeffs = CorrelatedRandomFieldIC::fourierCoeffs();

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


std::vector<std::vector<std::vector<std::vector<Real>>>>
CorrelatedRandomFieldIC::fourierCoeffs()
{
  int nodes = _Nnodes; //_mesh.n_nodes();
  std::array<Real, LIBMESH_DIM> L = {{_xmax - _xmin, _dim > 1 ? _ymax - _ymin : 0, _dim > 2 ? _zmax - _zmin : 0}};

  std::vector<std::vector<std::vector<std::vector<Real>>>> output = std::vector<std::vector<std::vector<std::vector<Real>>>>(2*nodes+1,std::vector<std::vector<std::vector<Real>>>(2*nodes+1,std::vector<std::vector<Real>>(2*nodes+1,std::vector<Real>(2,0))));

  Real rand1  = 0.0;
  Real rand2  = 0.0;

  Real nrand1  = 0.0;
  Real nrand2  = 0.0;

  Real randnorm = 0.0;
  Real amp = 0.0;

  RealVectorValue kv(0.0, 0.0, 0.0);
  Real k2 = 0.0;
  if (_dim == 1)
  {
    Real norm = std::sqrt(8.0*libMesh::pi/(_Lcorr*(_xmax-_xmin)));
    for (unsigned int i = 0; i<2*nodes+1; i++)
    {
      kv(0) = 2.0*libMesh::pi*double(i-nodes)/(_xmax-_xmin);
      k2 = kv(0)*kv(0); //+ kv(1)*kv(1) + kv(2)*kv(2);

      amp = norm/(k2 + Utility::pow<2>(1.0/_Lcorr));
      // Box-Muller method for gaussian random numbers

      rand1  = (rand()%1000000+1)/1000000.0;
      rand2  = (rand()%1000000+1)/1000000.0;

      //construct random numbers with unit variance
      //with unit variance for nrand1^2 + nrand2^2
      nrand1 = sqrt(-2.0*std::log(rand1))*cos(2.0*libMesh::pi*rand1)/sqrt(2.0);
      nrand2 = sqrt(-2.0*std::log(rand2))*sin(2.0*libMesh::pi*rand2)/sqrt(2.0);

      //enforce 0 mean
      if (i == nodes)
      {
        amp = 0.0;
      }
      output[i][0][0][0] = amp*nrand1;
      output[i][0][0][1] = amp*nrand2;
    }
  }
  else if (_dim == 2)
  {
    Real norm = std::sqrt(8.0*libMesh::pi/(_Lcorr*(_xmax-_xmin)*(_ymax-_ymin)));
    for (int i = 0; i<2*nodes+1; i++)
    {
      for (int j = 0; j<2*nodes+1; j++)
      {
        kv(0) = 2.0*libMesh::pi*double(i-nodes)/(_xmax-_xmin);
        kv(1) = 2.0*libMesh::pi*double(j-nodes)/(_ymax-_ymin);

        k2 = kv(0)*kv(0) + kv(1)*kv(1);

        amp = norm/(k2 + Utility::pow<2>(1.0/_Lcorr));
        // Box-Muller method for gaussian random numbers

        rand1  = (rand()%1000000+1)/1000000.0;
        rand2  = (rand()%1000000+1)/1000000.0;

        //construct random numbers with unit variance
        //with unit variance for nrand1^2 + nrand2^2
        nrand1 = sqrt(-2.0*std::log(rand1))*cos(2.0*libMesh::pi*rand1)/sqrt(2.0);
        nrand2 = sqrt(-2.0*std::log(rand2))*sin(2.0*libMesh::pi*rand2)/sqrt(2.0);

        //enforce 0 mean
        if (i == nodes && j == nodes)
        {
          amp = 0.0;
        }
        output[i][j][0][0] = amp*nrand1;
        output[i][j][0][1] = amp*nrand2;
      }
    }
  }
  else if (_dim == 3)
  {
    Real norm = std::sqrt(8.0*libMesh::pi/(_Lcorr*(_xmax-_xmin)*(_ymax-_ymin)*(_zmax-_zmin)));
    for (int i = 0; i<2*nodes+1; i++)
    {
      for (int j = 0; j<2*nodes+1; j++)
      {
        for (int k = 0; k<2*nodes+1; k++)
        {
          kv(0) = 2.0*libMesh::pi*double(i-nodes)/(_xmax-_xmin);
          kv(1) = 2.0*libMesh::pi*double(j-nodes)/(_ymax-_ymin);
          kv(2) = 2.0*libMesh::pi*double(k-nodes)/(_zmax-_zmin);

          k2 = kv(0)*kv(0) + kv(1)*kv(1) + kv(2)*kv(2);

          amp = norm/(k2 + Utility::pow<2>(1.0/_Lcorr));
          // Box-Muller method for gaussian random numbers

          rand1  = (rand()%1000000+1)/1000000.0;
          rand2  = (rand()%1000000+1)/1000000.0;

          //construct random numbers with unit variance
          //with unit variance for nrand1^2 + nrand2^2
          nrand1 = std::sqrt(-2.0*std::log(rand1))*cos(2.0*libMesh::pi*rand1)/std::sqrt(2.0);
          nrand2 = std::sqrt(-2.0*std::log(rand2))*sin(2.0*libMesh::pi*rand2)/std::sqrt(2.0);

          if (i == nodes && j == nodes && k == nodes)
          {
            amp = 0.0;
          }
          output[i][j][k][0] = amp*nrand1;
          output[i][j][k][1] = amp*nrand2;
        }
      }
    }
  }
  else
    Real norm = 0.0;
  return output;
}

Real
CorrelatedRandomFieldIC::evaluate(const Point & p, std::vector<std::vector<std::vector<std::vector<Real>>>> output)
{
  int nodes = _Nnodes;
  Real result = 0.0;
  RealVectorValue kv(0.0, 0.0, 0.0);
  Real phase = 0.0;

  if (_dim == 1)
  {
    for (unsigned int i = 0; i<2*nodes+1; i++)
    {
      kv(0) = 2.0*libMesh::pi*double(i-nodes)/(_xmax-_xmin);

      phase = kv(0)*p(0);
      //the real part
      result += output[i][0][0][0]*cos(phase);
      result += output[i][0][0][1]*sin(phase);
    }
  }
  else if (_dim == 2)
  {
    for (unsigned int i = 0; i<2*nodes+1; i++)
    {
      for (unsigned int j = 0;j<2*nodes+1;j++)
      {
        kv(0) = 2.0*libMesh::pi*double(i-nodes)/(_xmax-_xmin);
        kv(1) = 2.0*libMesh::pi*double(j-nodes)/(_ymax-_ymin);

        phase = kv(0)*p(0) + kv(1)*p(1) + kv(2)*p(2);
        //the real part
        result += output[i][j][0][0]*cos(phase);
        result += output[i][j][0][1]*sin(phase);
      }
    }
  }
  else if (_dim == 3)
  {
    for (int i = 0; i<2*nodes+1; i++)
    {
      for (int j = 0; j<2*nodes+1;j++)
      {
        for (int k = 0; k<2*nodes+1; k++)
        {
          kv(0) = 2.0*libMesh::pi*double(i-nodes)/(_xmax-_xmin);
          kv(1) = 2.0*libMesh::pi*double(j-nodes)/(_ymax-_ymin);
          kv(2) = 2.0*libMesh::pi*double(k-nodes)/(_zmax-_zmin);

          phase = kv(0)*p(0) + kv(1)*p(1) + kv(2)*p(2);
          //the real part

          result += output[i][j][k][0]*std::cos(phase);
          result += output[i][j][k][1]*std::sin(phase);
        }
      }
    }
  }
  else
    result = 0.0;
  return result;
}

Real
CorrelatedRandomFieldIC::value(const Point & p)
{
  return CorrelatedRandomFieldIC::evaluate(p, dcoeffs);
}
