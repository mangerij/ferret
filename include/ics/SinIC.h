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

#ifndef SINIC_H
#define SINIC_H

#include "Kernel.h"
#include "InitialCondition.h"
#include "InputParameters.h"

// Forward Declarations
class SinIC;

namespace libMesh { class Point; }

template<>
InputParameters validParams<SinIC>();


class SinIC:public InitialCondition
{
public:
  SinIC(const InputParameters & parameters);
  virtual Real value(const Point & p);
private:
  Real _amplitude;
  Real _wave_length_x,_wave_length_y,_wave_length_z;
  Real _phrase_x,_phrase_y,_phrase_z;
  Real _vertical_shift;
};

#endif //SinIC
