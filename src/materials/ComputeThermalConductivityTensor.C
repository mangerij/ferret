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

#include "ComputeThermalConductivityTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeThermalConductivityTensor);
InputParameters
ComputeThermalConductivityTensor::validParams()
{
  InputParameters params = ComputeRotatedThermalConductivityTensorBase::validParams();
  params.addClassDescription("Compute a ThermalConductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("k_ij", "ThermalConductivity tensor for material");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeThermalConductivityTensor::ComputeThermalConductivityTensor(
    const InputParameters & parameters) : ComputeRotatedThermalConductivityTensorBase(parameters)
    // _kij(this->template getParam<std::vector<Real> >("k_ij"), (RankTwoTensor::FillMethod)(int)this->template getParam<MooseEnum>("fill_method"))
{
  // _kij(getParam<std::vector<Real> >("k_ij"), (RankTwoTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"));
  _kij.fillFromInputVector(getParam<std::vector<Real>>("k_ij"));
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate thermal conductivity tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _kij.rotate(R);
}
void
ComputeThermalConductivityTensor::computeQpThermalConductivityTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _thC_tensor[_qp] = _kij;
}
