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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef MAGNETICANISOTROPYENERGY_H
#define MAGNETICANISOTROPYENERGY_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class MagneticAnisotropyEnergy;

template<>
InputParameters validParams<MagneticAnisotropyEnergy>();

//TODO: change the base class!
class MagneticAnisotropyEnergy : public ElementIntegralPostprocessor
{
public:
  MagneticAnisotropyEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const VariableValue& _mag_x;
  const VariableValue& _mag_y;
  const VariableValue& _mag_z;
  const Real _Ku;
  const Real _nx;
  const Real _ny;
  const Real _nz;
  const Real _len_scale;

};

#endif
