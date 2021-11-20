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

#ifndef TRANSFORMEDMICROFORCEELECTROSTRICTIVECOUPLINGENERGY_H
#define TRANFORMEDMICROFORCEELECTROSTRICTIVECOUPLINGENERGY_H

#include "AuxKernel.h"

class TransformedMicroforceElectrostrictiveCouplingEnergy : public AuxKernel
{
public:
  TransformedMicroforceElectrostrictiveCouplingEnergy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeValue();

private:
  const unsigned int _component;
  const unsigned int _u1_x_var;  // note that this is slightly different than other transformed 
  const unsigned int _u1_y_var;  // objects. We will set P = S^{-1} P' with the aux system but
  const unsigned int _u1_z_var;  // then u' = u' explicitly
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableGradient & _u1_x_grad;
  const VariableGradient & _u1_y_grad;
  const VariableGradient & _u1_z_grad;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _q11;
  const MaterialProperty<Real> & _q12;
  const MaterialProperty<Real> & _q44;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
};

#endif
