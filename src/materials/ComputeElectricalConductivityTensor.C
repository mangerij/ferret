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

#include "ComputeElectricalConductivityTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeElectricalConductivityTensor);
InputParameters
ComputeElectricalConductivityTensor::validParams()
{
  InputParameters params = ComputeRotatedElectricalConductivityTensorBase::validParams();
  params.addClassDescription("Compute a ElectricalConductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("g_ij", "ElectricalConductivity tensor for material");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeElectricalConductivityTensor::ComputeElectricalConductivityTensor(
    const InputParameters & parameters)
  : ComputeRotatedElectricalConductivityTensorBase(parameters)
// _gamma_ij(getParam<std::vector<Real> >("g_ij"),
// (RankTwoTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  _gij.fillFromInputVector(getParam<std::vector<Real>>("g_ij"));
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _gij.rotate(R);
}
void
ComputeElectricalConductivityTensor::computeQpElectricalConductivityTensor()
{
  ///Assign an electrical conductivity tensor at a given quad point. This will be reworked eventually for constant _qp.
  _ecC_tensor[_qp] = _gij;
}
