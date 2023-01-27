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

#include "ComputeRankFourLandauTensor.h"
#include "RotationTensor.h"
#include "RankFourTensor.h"

registerMooseObject("FerretApp", ComputeRankFourLandauTensor);

InputParameters ComputeRankFourLandauTensor::validParams()
{
  InputParameters params = ComputeRotatedRankFourLandauTensorBase::validParams();
  params.addClassDescription("Compute the quartic Landau coefficient tensor.");
  params.addRequiredParam<std::vector<Real> >("a_ijkl", "Landau tensor for material oriented [001]");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputeRankFourLandauTensor::ComputeRankFourLandauTensor(const InputParameters & parameters) :
    ComputeRotatedRankFourLandauTensorBase(parameters),
    _alpha_ijkl(getParam<std::vector<Real> >("a_ijkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles);
  _alpha_ijkl.rotate(R);
}

void
ComputeRankFourLandauTensor::computeQpRankFourLandauTensor()
{
  _rank_four_landau_tensor[_qp] = _alpha_ijkl;
}
