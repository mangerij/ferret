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
   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret
**/

#ifndef TPELECTROSTRICTIVECOUPLINGPDERIVATIVE_H
#define TPELECTROSTRICTIVECOUPLINGPDERIVATIVE_H

#include "Kernel.h"

class TPElectrostrictiveCouplingPDerivative: public Kernel
{
public:
  TPElectrostrictiveCouplingPDerivative(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _C11;
  const MaterialProperty<Real> & _C12;
  const MaterialProperty<Real> & _C13;
  const MaterialProperty<Real> & _C33;
  const MaterialProperty<Real> & _C44;
  const MaterialProperty<Real> & _C66;
  const MaterialProperty<Real> & _Q11;
  const MaterialProperty<Real> & _Q12;
  const MaterialProperty<Real> & _Q13;
  const MaterialProperty<Real> & _Q33;
  const MaterialProperty<Real> & _Q44;
  const MaterialProperty<Real> & _Q66;
  const std::string _base_name;
  const MaterialProperty<RankTwoTensor> & _strain;
};
#endif //TPELECTROSTRICTIVECOUPLINGPDERIVATIVE_H
