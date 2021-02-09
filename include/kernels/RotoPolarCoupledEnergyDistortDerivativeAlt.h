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

#ifndef ROTOPOLARCOUPLEDENERGYDISTORTDERIVATIVEALT_H
#define ROTOPOLARCOUPLEDENERGYDISTORTDERIVATIVEALT_H

#include "Kernel.h"

class RotoPolarCoupledEnergyDistortDerivativeAlt;

template<>
InputParameters validParams<RotoPolarCoupledEnergyDistortDerivativeAlt>();

class RotoPolarCoupledEnergyDistortDerivativeAlt: public Kernel
{
public:

  RotoPolarCoupledEnergyDistortDerivativeAlt(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _antiferrodis_A_x_var;
  const unsigned int _antiferrodis_A_y_var;
  const unsigned int _antiferrodis_A_z_var;
  const VariableValue & _antiferrodis_A_x;
  const VariableValue & _antiferrodis_A_y;
  const VariableValue & _antiferrodis_A_z;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _t1111;
  const MaterialProperty<Real> & _t1122;
  const MaterialProperty<Real> & _t1212;
  const MaterialProperty<Real> & _t42111111;
  const MaterialProperty<Real> & _t24111111;
  const MaterialProperty<Real> & _t42111122;
  const MaterialProperty<Real> & _t24112222;
  const MaterialProperty<Real> & _t42112233;
  const MaterialProperty<Real> & _t24112233;
  const MaterialProperty<Real> & _t42112211;
  const MaterialProperty<Real> & _t24111122;
  const MaterialProperty<Real> & _t42111212;
  const MaterialProperty<Real> & _t42123312;
  const MaterialProperty<Real> & _t24121112;
  const MaterialProperty<Real> & _t24121233;
  const MaterialProperty<Real> & _t6211111111;
  const MaterialProperty<Real> & _t2611111111;
  const MaterialProperty<Real> & _t6211111122;
  const MaterialProperty<Real> & _t2611222222;
  const MaterialProperty<Real> & _t4411111111;
  const MaterialProperty<Real> & _t4411112222;
};
#endif //ROTOPOLARCOUPLEDENERGYDISTORTDERIVATIVEALT_H
