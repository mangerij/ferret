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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "MagnetostrictiveCouplingHeff.h"
#include "libmesh/utility.h"

class MagnetostrictiveCouplingHeff;

registerMooseObject("FerretApp", MagnetostrictiveCouplingHeff);

template<>
InputParameters validParams<MagnetostrictiveCouplingHeff>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to the magnetoelectric effective field. Note for cubic magnets only.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization");
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elastic displacement");
  return params;
}

MagnetostrictiveCouplingHeff::MagnetostrictiveCouplingHeff(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _disp_x(coupledValue("disp_x")),
   _disp_y(coupledValue("disp_y")),
   _disp_z(coupledValue("disp_z")),
   _alpha(getMaterialProperty<Real>("alpha")),
   _g0(getMaterialProperty<Real>("g0")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _L11(getMaterialProperty<Real>("L11")),
   _L12(getMaterialProperty<Real>("L12")),
   _L44(getMaterialProperty<Real>("L44"))
{
}

Real
MagnetostrictiveCouplingHeff::computeQpResidual()
{
  if (_component == 0)
  {
    return 0.0;
  }
  else if (_component == 1)
  {
    return 0.0;
  }
  else if (_component == 2)
  {
    return 0.0;
  }
  else
    return 0.0;
}

Real
MagnetostrictiveCouplingHeff::computeQpJacobian()
{
  return 0.0;
}

Real
MagnetostrictiveCouplingHeff::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
