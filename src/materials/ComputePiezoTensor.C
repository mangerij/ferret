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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/
#include "ComputeElasticityTensor.h"
#include "ComputeRotatedElasticityTensorBase.h"
#include "ComputePiezoTensor.h"
#include "RotationTensor.h"
#include "RankThreeTensor.h"
#include "PiezostrictiveTensorTools.h"

template<>
InputParameters validParams<ComputePiezoTensor>()
{
  InputParameters params = validParams<ComputeRotatedPiezoTensorBase>();
  params.addClassDescription("Compute the converse piezoelectric tensor.");
  params.addParam<bool>("compute_piezostrictive_coeff", false, "compute the piezostrictive coefficients C_ijmn * d_kij");
  params.addRequiredParam<std::vector<Real> >("d_kpq", "piezoelectric tensor for material");
  params.addRequiredParam<std::vector<Real> >("d_pqkT", "piezoelectric tensor for material");
  params.addRequiredParam<std::vector<Real> >("C_ijkl", "elastic stiffness tensor for material");
  params.addParam<MooseEnum>("fill_method", RankThreeTensor::fillMethodEnum() = "general", "The fill method");
  params.addParam<MooseEnum>("fill_method", RankThreeTensor::fillMethodEnum() = "general", "The fill method");
  params.addParam<MooseEnum>("fill_method2", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputePiezoTensor::ComputePiezoTensor(const InputParameters & parameters) :
    ComputeRotatedPiezoTensorBase(parameters),
    _dkpq(getParam<std::vector<Real> >("d_kpq"), (RankThreeTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _dpqkT(getParam<std::vector<Real> >("d_pqkT"), (RankThreeTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _Cijkl(getParam<std::vector<Real> >("C_ijkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method2"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  // _dmkl.rotate(R);
  ///contractions using namespace method
  if (_compute_piezostrictive_coeff == true)
    _Dkij = PiezostrictiveTensorTools::computeProduct(_Cijkl, _dkpq);
    Moose::out << " \n     ";
    Moose::out << "_Dkij(0,0,0) "; std::cout << _Dkij(0,0,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,0,1) "; std::cout << _Dkij(0,0,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,0,2) "; std::cout << _Dkij(0,0,2); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,1,0) "; std::cout << _Dkij(0,1,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,1,1) "; std::cout << _Dkij(0,1,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,1,2) "; std::cout << _Dkij(0,1,2); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,2,0) "; std::cout << _Dkij(0,2,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,2,1) "; std::cout << _Dkij(0,2,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(0,2,2) "; std::cout << _Dkij(0,2,2); Moose::out << " \n     ";
    Moose::out << " \n     ";
    Moose::out << "_Dkij(1,0,0) "; std::cout << _Dkij(1,0,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,0,1) "; std::cout << _Dkij(1,0,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,0,2) "; std::cout << _Dkij(1,0,2); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,1,0) "; std::cout << _Dkij(1,1,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,1,1) "; std::cout << _Dkij(1,1,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,1,2) "; std::cout << _Dkij(1,1,2); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,2,0) "; std::cout << _Dkij(1,2,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,2,1) "; std::cout << _Dkij(1,2,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(1,2,2) "; std::cout << _Dkij(1,2,2); Moose::out << " \n     ";
    Moose::out << " \n     ";
    Moose::out << "_Dkij(2,0,0) "; std::cout << _Dkij(2,0,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,0,1) "; std::cout << _Dkij(2,0,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,0,2) "; std::cout << _Dkij(2,0,2); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,1,0) "; std::cout << _Dkij(2,1,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,1,1) "; std::cout << _Dkij(2,1,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,1,2) "; std::cout << _Dkij(2,1,2); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,2,0) "; std::cout << _Dkij(2,2,0); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,2,1) "; std::cout << _Dkij(2,2,1); Moose::out << " \n     ";
    Moose::out << "_Dkij(2,2,2) "; std::cout << _Dkij(2,2,2); Moose::out << " \n     ";
    _DijkT = PiezostrictiveTensorTools::computePiezoTransposeProduct(_Cijkl, _dpqkT);
    Moose::out << " \n     ";
    Moose::out << "_DijkT(0,0,0)"; std::cout << _DijkT(0,0,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,0,1)"; std::cout << _DijkT(0,0,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,0,2)"; std::cout << _DijkT(0,0,2); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,1,0)"; std::cout << _DijkT(0,1,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,1,1)"; std::cout << _DijkT(0,1,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,1,2)"; std::cout << _DijkT(0,1,2); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,2,0)"; std::cout << _DijkT(0,2,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,2,1)"; std::cout << _DijkT(0,2,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(0,2,2)"; std::cout << _DijkT(0,2,2); Moose::out << " \n     ";
    Moose::out << " \n     ";
    Moose::out << "_DijkT(1,0,0)"; std::cout << _DijkT(1,0,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,0,1)"; std::cout << _DijkT(1,0,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,0,2)"; std::cout << _DijkT(1,0,2); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,1,0)"; std::cout << _DijkT(1,1,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,1,1)"; std::cout << _DijkT(1,1,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,1,2)"; std::cout << _DijkT(1,1,2); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,2,0)"; std::cout << _DijkT(1,2,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,2,1)"; std::cout << _DijkT(1,2,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(1,2,2)"; std::cout << _DijkT(1,2,2); Moose::out << " \n     ";
    Moose::out << " \n     ";
    Moose::out << "_DijkT(2,0,0)"; std::cout << _DijkT(2,0,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,0,1)"; std::cout << _DijkT(2,0,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,0,2)"; std::cout << _DijkT(2,0,2); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,1,0)"; std::cout << _DijkT(2,1,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,1,1)"; std::cout << _DijkT(2,1,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,1,2)"; std::cout << _DijkT(2,1,2); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,2,0)"; std::cout << _DijkT(2,2,0); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,2,1)"; std::cout << _DijkT(2,2,1); Moose::out << " \n     ";
    Moose::out << "_DijkT(2,2,2)"; std::cout << _DijkT(2,2,2); Moose::out << " \n     ";
    Moose::out << " \n     ";
}

void
ComputePiezoTensor::computeQpPiezoTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _piezo_tensor[_qp] = _dkpq;
  if (_compute_piezostrictive_coeff == true)
    _piezostrictive_tensor[_qp] = _Dkij;
    _piezostrictive_tensor_i[_qp] = _DijkT;

}
