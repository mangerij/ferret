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

#ifndef MAGNETICEXCHANGEENERGYDENSITYUS_H
#define MAGNETICEXCHANGEENERGYDENSITYUS_H

#include "AuxKernel.h"

class MagneticExchangeEnergyDensityUS : public AuxKernel
{
public:
  MagneticExchangeEnergyDensityUS(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeValue();
  const VariableValue & _azimuthal_ph;
  const VariableValue & _polar_th;
  const VariableGradient & _polar_th_grad;
  const VariableGradient & _azimuthal_ph_grad;
  const MaterialProperty<Real> & _Ae;
  const MaterialProperty<Real> & _Ms;
  const Real _energy_scale;
};

#endif // MAGNETICEXCHANGEENERGYDENSITYUS_H
