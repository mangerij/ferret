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

#ifndef INHOMOGENEOUSBULKENERGY_H
#define INHOMOGENEOUSBULKENERGY_H

#include "ElementIntegralPostprocessor.h"

class InhomogeneousBulkEnergy : public ElementIntegralPostprocessor
{
public:
  InhomogeneousBulkEnergy(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual Real computeQpIntegral();
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableValue & _x;
  const Real _alpha01, _alpha011, _alpha0111, _alpha012, _alpha0112, _alpha0123, _alpha1111, _alpha1112, _alpha1122, _alpha1123, _b1, _b2, _b3, _T;
  const Real _len_scale;

};

#endif
