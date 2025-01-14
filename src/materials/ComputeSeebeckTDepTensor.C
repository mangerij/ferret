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

#include "ComputeSeebeckTDepTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeSeebeckTDepTensor);
InputParameters
ComputeSeebeckTDepTensor::validParams()
{
  InputParameters params = ComputeRotatedSeebeckTensorBase::validParams();
  params.addClassDescription("Compute a Seebeck tensor.");
  params.addRequiredParam<std::vector<Real>>("asb_ij", "Seebeck tensor for material");
  params.addRequiredParam<std::vector<Real>>("bsb_ij", "Seebeck tensor for material");
  params.addRequiredParam<std::vector<Real>>("csb_ij", "Seebeck tensor for material");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeSeebeckTDepTensor::ComputeSeebeckTDepTensor(const InputParameters & parameters)
  : ComputeRotatedSeebeckTensorBase(parameters),
   _T(coupledValue("T"))
{
  _asbij.fillFromInputVector(getParam<std::vector<Real>>("asb_ij"));
  _bsbij.fillFromInputVector(getParam<std::vector<Real>>("bsb_ij"));
  _csbij.fillFromInputVector(getParam<std::vector<Real>>("csb_ij"));
  RotationTensor R(_Euler_angles);

  _asbij.rotate(R);
  _bsbij.rotate(R);
  _csbij.rotate(R);
}

void
ComputeSeebeckTDepTensor::computeQpSeebeckTensor()
{
  ///Assign a Seebeck tensor at a given quad point. This will be reworked eventually for constant _qp.
  _sbC_tensor[_qp] = _asbij+_bsbij*_T[_qp]+_csbij*_T[_qp]*_T[_qp];
}
