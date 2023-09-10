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

#include "ComputePiezostrictiveTensor.h"
#include "RotationTensor.h"
#include "RankThreeTensor.h"

registerMooseObject("FerretApp", ComputePiezostrictiveTensor);

InputParameters ComputePiezostrictiveTensor::validParams()
{
  InputParameters params = ComputeRotatedPiezostrictiveTensorBase::validParams();
  params.addClassDescription("Compute a piezostrictive tensor.");
  params.addRequiredParam<std::vector<Real> >("e_ijk", "piezostrictive tensor for material");
  params.addParam<MooseEnum>("fill_method", RankThreeTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}

ComputePiezostrictiveTensor::ComputePiezostrictiveTensor(const InputParameters & parameters) :
    ComputeRotatedPiezostrictiveTensorBase(parameters),
    _eijk(getParam<std::vector<Real> >("e_ijk"), (RankThreeTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _eijk.rotate(R);
}

void
ComputePiezostrictiveTensor::computeQpPiezostrictiveTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _piezostrictive_tensor[_qp] = _eijk;
}
