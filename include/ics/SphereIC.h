/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#ifndef SPHEREIC_H
#define SPHEREIC_H

#include "Kernel.h"
#include "InitialCondition.h"
#include "InputParameters.h"
// Forward Declarations
class SphereIC;
class Function;

namespace libMesh { class Point; }

template<>
InputParameters validParams<SphereIC>();


class SphereIC:public InitialCondition
{
public:
  SphereIC(const InputParameters & parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);
protected:
  Function & _radial_func;
  Function & _polar_func;
  Function & _azimuthal_func;
  unsigned int _index;
};

#endif //SphereIC
