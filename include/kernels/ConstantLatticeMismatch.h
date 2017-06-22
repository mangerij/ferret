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

#ifndef CONSTANTLATTICEMISMATCH_H
#define CONSTANTLATTICEMISMATCH_H

#include "Kernel.h"

class ConstantLatticeMismatch;

template<>
InputParameters validParams<ConstantLatticeMismatch>();

class ConstantLatticeMismatch: public Kernel
{
public:

  ConstantLatticeMismatch(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const VariableValue & _Qxx;
  const VariableValue & _Qxy;
  const VariableValue & _Qxz;
  const VariableValue & _Qyx;
  const VariableValue & _Qyy;
  const VariableValue & _Qyz;
  const VariableValue & _Qzx;
  const VariableValue & _Qzy;
  const VariableValue & _Qzz;
  const unsigned int _disp_x_var;
  const unsigned int _disp_y_var;
  const unsigned int _disp_z_var;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const unsigned int _component;
  const unsigned int _deriv_component;
};
#endif
