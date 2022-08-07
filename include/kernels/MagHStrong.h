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

#ifndef MAGHSTRONG_H
#define MAGHSTRONG_H

#include "Kernel.h"

class MagHStrong: public Kernel
{
public:

  MagHStrong(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _azimuthal_ph_var;
  const unsigned int _polar_th_var;
  const VariableValue & _azimuthal_ph;
  const VariableValue & _polar_th;
  const MaterialProperty<Real> & _mu0;
  const MaterialProperty<Real> & _Ms;

};
#endif //MAGHSTRONG_H
