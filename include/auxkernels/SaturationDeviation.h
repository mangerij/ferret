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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#ifndef SATURATIONDEVIATION_H
#define SATURATIONDEVIATION_H

#include "AuxKernel.h"

class SaturationDeviation;

template<>
InputParameters validParams<SaturationDeviation>();

class SaturationDeviation: public AuxKernel
{
public:
  SaturationDeviation(const InputParameters & parameters);

  virtual ~SaturationDeviation() {}

protected:
    virtual Real computeValue();

private:
  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  const Real _Ms;
};

#endif // SATURATIONDEVIATION_H