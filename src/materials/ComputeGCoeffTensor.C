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

#include "ComputeGCoeffTensor.h"
#include "RotationTensor.h"
#include "RankFourTensor.h"

template<>
InputParameters validParams<ComputeGCoeffTensor>()
{
  InputParameters params = validParams<ComputeRotatedGCoeffTensorBase>();
  params.addClassDescription("Compute a polar-optic (g) tensor.");
  params.addRequiredParam<std::vector<Real> >("g_ijkl", "polar-optic tensor for material");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputeGCoeffTensor::ComputeGCoeffTensor(const InputParameters & parameters) :
    ComputeRotatedGCoeffTensorBase(parameters),
    _gijkl(getParam<std::vector<Real> >("g_ijkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate polar-optic tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _gijkl.rotate(R);
}

void
ComputeGCoeffTensor::computeQpGCoeffTensor()
{
  ///Assign a polar-optic tensor at a given quad point. This will be reworked eventually for constant _qp.
  _gcoefficient_tensor[_qp] = _gijkl;
}

