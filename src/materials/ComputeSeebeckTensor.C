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

#include "ComputeSeebeckTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeSeebeckTensor);
InputParameters
ComputeSeebeckTensor::validParams()
{
  InputParameters params = ComputeRotatedSeebeckTensorBase::validParams();
  params.addClassDescription("Compute a Seebeck tensor.");
  params.addRequiredParam<std::vector<Real>>("a_ij", "Seebeck tensor for material");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeSeebeckTensor::ComputeSeebeckTensor(const InputParameters & parameters)
  : ComputeRotatedSeebeckTensorBase(parameters)
// _alpha_ij(getParam<std::vector<Real> >("a_ij"),
// (RankTwoTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  _aij.fillFromInputVector(getParam<std::vector<Real>>("a_ij"));
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate seebeck tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _aij.rotate(R);
}
void
ComputeSeebeckTensor::computeQpSeebeckTensor()
{
  ///Assign a seebeck tensor at a given quad point. This will be reworked eventually for constant _qp.
  _sbC_tensor[_qp] = _aij;
}
