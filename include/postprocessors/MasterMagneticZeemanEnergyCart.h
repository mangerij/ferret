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

#ifndef MASTERMAGNETICZEEMANENERGYCART_H
#define MASTERMAGNETICZEEMANENERGYCART_H

#include "ElementIntegralPostprocessor.h"

class MasterMagneticZeemanEnergyCart : public ElementIntegralPostprocessor
{
public:
  MasterMagneticZeemanEnergyCart(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpIntegral();
  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  
  const VariableValue &  _Hext_x;
  const VariableValue &  _Hext_y;
  const VariableValue &  _Hext_z;
  const Real & _energy_scale;
  const MaterialProperty<Real> & _Ms;
  const Real & _Hscale;
  const Real & _mu0;

};

#endif
