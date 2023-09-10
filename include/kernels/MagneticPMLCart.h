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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#ifndef MAGNETICPMLCART_H
#define MAGNETICPMLCART_H

#include "Kernel.h"

class MagneticPMLCart: public Kernel
{
public:

  MagneticPMLCart(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const unsigned int _component;
  const MaterialProperty<Real> & _permittivity;
  /*   const Real &   _deltasxminus;
   const Real &   _deltapxminus;
   const Real &   _deltawxminus;
   const Real &   _x0pmlminus;*/
   const Real &   _deltasyminus;
   const Real &   _deltapyminus;
   const Real &   _deltawyminus;
   const Real &   _y0pmlminus;
  /*   const Real &   _deltaszminus;
   const Real &   _deltapzminus;
   const Real &   _deltawzminus;
   const Real &   _z0pmlminus;
   const Real &   _deltasxplus;
   const Real &   _deltapxplus;
   const Real &   _deltawxplus;
   const Real &   _x0pmlplus;
   const Real &   _deltasyplus;
   const Real &   _deltapyplus;
   const Real &   _deltawyplus;
   const Real &   _y0pmlplus;
   const Real &   _deltaszplus;
   const Real &   _deltapzplus;
   const Real &   _deltawzplus;
   const Real &   _z0pmlplus;*/
  

};
#endif
