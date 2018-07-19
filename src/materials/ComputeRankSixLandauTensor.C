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

#include "ComputeRankSixLandauTensor.h"
#include "RotationTensor.h"
#include "RankSixTensor.h"

registerMooseObject("FerretApp", ComputeRankSixLandauTensor);

template<>
InputParameters validParams<ComputeRankSixLandauTensor>()
{
  InputParameters params = validParams<ComputeRotatedRankSixLandauTensorBase>();
  params.addClassDescription("Compute the sextic Landau coefficient tensor.");
  params.addRequiredParam<std::vector<Real> >("a_ijklmn", "Landau tensor for material oriented [001]");
  params.addParam<MooseEnum>("fill_method", RankSixTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}

ComputeRankSixLandauTensor::ComputeRankSixLandauTensor(const InputParameters & parameters) :
    ComputeRotatedRankSixLandauTensorBase(parameters),
    _alpha_ijklmn(getParam<std::vector<Real> >("a_ijklmn"), (RankSixTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles);
  _alpha_ijklmn.rotate(R);
}

void
ComputeRankSixLandauTensor::computeQpRankSixLandauTensor()
{
  _rank_six_landau_tensor[_qp] = _alpha_ijklmn;
}
