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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ComputeThermalConductivityTDepTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeThermalConductivityTDepTensor);
InputParameters
ComputeThermalConductivityTDepTensor::validParams()
{
  InputParameters params = ComputeRotatedThermalConductivityTensorBase::validParams();
  params.addClassDescription("Compute a ThermalConductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("ak_ij", "ThermalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("bk_ij", "ThermalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("ck_ij", "ThermalConductivity tensor for material");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeThermalConductivityTDepTensor::ComputeThermalConductivityTDepTensor(
    const InputParameters & parameters) : ComputeRotatedThermalConductivityTensorBase(parameters),
   _T(coupledValue("T"))
{
  _akij.fillFromInputVector(getParam<std::vector<Real>>("ak_ij"));
  _bkij.fillFromInputVector(getParam<std::vector<Real>>("bk_ij"));
  _ckij.fillFromInputVector(getParam<std::vector<Real>>("ck_ij"));

  RotationTensor R(_Euler_angles);
  _akij.rotate(R);
  _bkij.rotate(R);
  _ckij.rotate(R);
}
void
ComputeThermalConductivityTDepTensor::computeQpThermalConductivityTensor()
{
  _thC_tensor[_qp] = _akij + _bkij*_T[_qp] + _ckij*_T[_qp]*_T[_qp];
}
