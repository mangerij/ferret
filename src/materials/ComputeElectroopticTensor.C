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

**/

#include "ComputeElectroopticTensor.h"
#include "RotationTensor.h"
#include "RankThreeTensor.h"

registerMooseObject("FerretApp", ComputeElectroopticTensor);

template<>
InputParameters validParams<ComputeElectroopticTensor>()
{
  InputParameters params = validParams<ComputeRotatedElectroopticTensorBase>();
  params.addClassDescription("Compute an electrooptic tensor.");
  params.addRequiredParam<std::vector<Real> >("r_ijk", "electrooptic tensor for material");
  params.addParam<MooseEnum>("fill_method", RankThreeTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}

ComputeElectroopticTensor::ComputeElectroopticTensor(const InputParameters & parameters) :
    ComputeRotatedElectroopticTensorBase(parameters),
    _rijk(getParam<std::vector<Real> >("r_ijk"), (RankThreeTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _rijk.rotate(R);
}

void
ComputeElectroopticTensor::computeQpElectroopticTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _electrooptic_tensor[_qp] = _rijk;
}

