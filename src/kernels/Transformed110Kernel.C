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

#include "Transformed110Kernel.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", Transformed110Kernel);

template<>
InputParameters validParams<Transformed110Kernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the transformed residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2.0 for z)");
  params.addRequiredCoupledVar("order_param_x", "The x component of the transformed order parameter");
  params.addRequiredCoupledVar("order_param_y", "The y component of the transformed order parameter");
  params.addCoupledVar("order_param_z", 0.0, "The z component of the transformed order parameter");
  params.addRequiredCoupledVar("fb_x", "The x component of the untransformed force");
  params.addRequiredCoupledVar("fb_y", "The y component of the untransformed force");
  params.addRequiredCoupledVar("fb_z", "The z component of the untransformed force");
  params.addRequiredCoupledVar("Jb_xx", "The xx component of the untransformed jacobian");
  params.addRequiredCoupledVar("Jb_yy", "The yy component of the untransformed jacobian");
  params.addRequiredCoupledVar("Jb_zz", "The zz component of the untransformed jacobian");
  params.addRequiredCoupledVar("Jb_xy", "The xy component of the untransformed jacobian");
  params.addRequiredCoupledVar("Jb_yz", "The yz component of the untransformed jacobian");
  params.addRequiredCoupledVar("Jb_xz", "The xz component of the untransformed jacobian");
  return params;
}

Transformed110Kernel::Transformed110Kernel(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _order_param_x_var(coupled("order_param_x")),
   _order_param_y_var(coupled("order_param_y")),
   _order_param_z_var(coupled("order_param_z")),
   _fb_x(coupledValue("fb_x")),
   _fb_y(coupledValue("fb_y")),
   _fb_z(coupledValue("fb_z")),
   _Jb_xx(coupledValue("Jb_xx")),
   _Jb_yy(coupledValue("Jb_yy")),
   _Jb_zz(coupledValue("Jb_zz")),
   _Jb_xy(coupledValue("Jb_xy")),
   _Jb_yz(coupledValue("Jb_yz")),
   _Jb_xz(coupledValue("Jb_xz"))
{
}

Real
Transformed110Kernel::computeQpResidual()

//
// TODO: Note that there is no reason this needs to be hardcoded, but this will be the first step. 
//       in general, this procedure should work for any transformation
//
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (0.0);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (0.0);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (0.0);
  }
  else
    return 0.0;
}

Real
Transformed110Kernel::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
  }
  else
    return 0.0;
}


Real
Transformed110Kernel::computeQpOffDiagJacobian(unsigned int jvar)
{ 
  if (_component == 0)
  {
    if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.0);
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.0);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
    }
    else if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
