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

#ifndef BULKENERGYCOUPLEDT_H
#define BULKENERGYCOUPLEDT_H

#include "ElementIntegralPostprocessor.h"

class BulkEnergyCoupledT;

template<>
InputParameters validParams<BulkEnergyCoupledT>();

class BulkEnergyCoupledT : public ElementIntegralPostprocessor
{
public:
  BulkEnergyCoupledT(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const VariableValue & _temperature;
  const Real _alpha0, _alpha11, _alpha12, _alpha111, _alpha112, _alpha123, _Tc;
  const Real _len_scale;

};

#endif
