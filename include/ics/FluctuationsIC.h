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

#ifndef FLUCTUATIONSIC_H
#define FLUCTUATIONSIC_H

#include "InitialCondition.h"

// Forward Declarations
class FluctuationsIC;
namespace libMesh { class Point; }

template<>
InputParameters validParams<FluctuationsIC>();

class FluctuationsIC : public InitialCondition
{
public:
  FluctuationsIC(const InputParameters & parameters);
  virtual Real value(const Point & p);

protected:

private:
  Real _epsilon;
  Real _base_value;

  Point _q1;
  Point _q2;
  Point _q3;
  Point _q4;

  Real _h;

};

#endif //FLUCTUATIONSIC_H
