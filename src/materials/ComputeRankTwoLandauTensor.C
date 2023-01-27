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

#include "ComputeRankTwoLandauTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeRankTwoLandauTensor);

InputParameters ComputeRankTwoLandauTensor::validParams()
{
  InputParameters params = ComputeRotatedRankTwoLandauTensorBase::validParams();
  params.addClassDescription("Compute the quadratic Landau coefficient tensor.");
  params.addRequiredParam<std::vector<Real> >("a_ij", "Landau tensor for material oriented [001]");
  //params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "symmetric6", "The fill method");
  return params;
}

ComputeRankTwoLandauTensor::ComputeRankTwoLandauTensor(const InputParameters & parameters) :
    ComputeRotatedRankTwoLandauTensorBase(parameters)
{
  _alpha_ij.fillFromInputVector(getParam<std::vector<Real>>("a_ij"));
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles);
  _alpha_ij.rotate(R);
}

void
ComputeRankTwoLandauTensor::computeQpRankTwoLandauTensor()
{
  _rank_two_landau_tensor[_qp] = _alpha_ij;
}
