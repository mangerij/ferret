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

#include "TensorPressureAux.h"

registerMooseObject("FerretApp", TensorPressureAux);

InputParameters TensorPressureAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the value of the hydrostatic stress (which is 1/3 the minus of the stress tensor trace).");
  return params;
}

TensorPressureAux::TensorPressureAux(const InputParameters & parameters) :
  AuxKernel(parameters ),
   _stress( getMaterialProperty<RankTwoTensor>("stress") )
{}

Real
TensorPressureAux::computeValue()
{
    return -0.33333333333*_stress[_qp].trace();

}


