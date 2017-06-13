/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "ComputeDeltaIndicatrixElectro.h"
#include "RankThreeTensor.h"

template<>
InputParameters validParams<ComputeDeltaIndicatrixElectro>()
{
  InputParameters params = validParams<ComputeElectroopticTensorBase>();
  params.addClassDescription("Compute the adjustments to the indicatrix (beta tensor).");
  return params;
}

ComputeDeltaIndicatrixElectro::ComputeDeltaIndicatrixElectro(const InputParameters & parameters) :
    ComputeElectroopticTensorBase(parameters),
    _permittivity(getMaterialProperty<Real>("permittivity")),
    _len_scale(getMaterialProperty<Real>("len_scale")),
    _electrooptic_tensor(getMaterialProperty<RankFourTensor>("electrooptic_tensor"))
{
}

void
ComputeDeltaIndicatrixElectro::computeQpDeltaIndicatrix()
{
  Real sum = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      for (unsigned int k = 0; k < 3; ++k)
        {
          sum += _electrooptic_tensor[_qp](i, j, k) * (_permittivity * _grad_potential_int[_qp](k) * _len_scale);
        }
    _delta_indicatrix[_qp](i) = sum;
    }
    //Moose::out << "\n b"; std::cout << a; Moose::out << " = "; std::cout << _delta_beta_tensor[_qp](0, a);
}
