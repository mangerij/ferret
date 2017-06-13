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

#ifndef SURFACECHARGEAUX_H
#define SURFACECHARGEAUX_H

#include "AuxKernel.h"


//Forward declarations
class SurfaceChargeAux;

template<>
InputParameters validParams<SurfaceChargeAux>();

class SurfaceChargeAux : public AuxKernel
{
public:
  SurfaceChargeAux(const InputParameters & parameters);

  virtual ~SurfaceChargeAux() {}

protected:
  virtual Real computeValue();
  const MooseArray<Point> & _normals;

private:
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
};

#endif
