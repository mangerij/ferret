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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ComputeElectricalConductivityTDepTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeElectricalConductivityTDepTensor);
InputParameters
ComputeElectricalConductivityTDepTensor::validParams()
{
  InputParameters params = ComputeRotatedElectricalConductivityTensorBase::validParams();
  params.addClassDescription("Store a temperature dependent electrical conductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("asg_ij", "ElectricalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("bsg_ij", "ElectricalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("csg_ij", "ElectricalConductivity tensor for material");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeElectricalConductivityTDepTensor::ComputeElectricalConductivityTDepTensor(
    const InputParameters & parameters)
  : ComputeRotatedElectricalConductivityTensorBase(parameters),
   _T(coupledValue("T"))
{
  _asgij.fillFromInputVector(getParam<std::vector<Real>>("asg_ij"));
  _bsgij.fillFromInputVector(getParam<std::vector<Real>>("bsg_ij"));
  _csgij.fillFromInputVector(getParam<std::vector<Real>>("csg_ij"));

  RotationTensor R(_Euler_angles);

  _asgij.rotate(R);
  _bsgij.rotate(R);
  _csgij.rotate(R);
}

void
ComputeElectricalConductivityTDepTensor::computeQpElectricalConductivityTensor()
{
  ///Assign an electrical conductivity tensor at a given quad point. This will be reworked eventually for constant _qp.
  _ecC_tensor[_qp] = _asgij+_bsgij*_T[_qp] + _csgij*_T[_qp]*_T[_qp];
}
