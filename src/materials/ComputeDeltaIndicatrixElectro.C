/**************************************************************************
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

***************************************************************************/

#include "ComputeDeltaIndicatrixElectro.h"
#include "RankThreeTensor.h"
#include "ComputeElectroopticTensor.h"

template<>
InputParameters validParams<ComputeDeltaIndicatrixElectro>()
{
  InputParameters params = validParams<ComputeDeltaIndicatrixElectroBase>();
  params.addClassDescription("Compute the adjustments to the indicatrix (beta tensor).");
  params.addRequiredCoupledVar("potential_E_int", "The electrostatic potential");
  return params;
}

ComputeDeltaIndicatrixElectro::ComputeDeltaIndicatrixElectro(const InputParameters & parameters) :
    ComputeDeltaIndicatrixElectroBase(parameters),
    _electrooptic_tensor(getMaterialProperty<RankThreeTensor>("electrooptic_tensor")),
    _potential_E_int_grad(coupledGradient("potential_E_int"))
{
}

void
ComputeDeltaIndicatrixElectro::computeQpDeltaIndicatrixElectro()
{
  Real sum = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      for (unsigned int k = 0; k < 3; ++k)
        {
          sum += _electrooptic_tensor[_qp](i, j, k) * _potential_E_int_grad[_qp](k) ;
        }
    _delta_indicatrix_electro[_qp](i, j) = sum;
    }
}
