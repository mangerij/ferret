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

#ifndef MAGNETICANISOTROPYENERGY_H
#define MAGNETICANISOTROPYENERGY_H

#include "ElementIntegralPostprocessor.h"

class MagneticAnisotropyEnergy : public ElementIntegralPostprocessor
{
public:
  MagneticAnisotropyEnergy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpIntegral();
  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  const MaterialProperty<Real> & _K1;
  const MaterialProperty<Real> & _K2;
  const MaterialProperty<Real> & _nx;
  const MaterialProperty<Real> & _ny;
  const MaterialProperty<Real> & _nz;
//  const MaterialProperty<Real> & _Ms;
//  const MaterialProperty<Real> & _mu0;
};

#endif
