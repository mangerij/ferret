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

   You should have received a co_antiferrodis_A_y[_i] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ScalarRotostrictiveEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ScalarRotostrictiveEnergy);

InputParameters
ScalarRotostrictiveEnergy::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addCoupledVar("antiferrodis_A_x", "variable antiferrodis_A_x coupled into this kernel");
  params.addCoupledVar("antiferrodis_A_y", "variable antiferrodis_A_y coupled into this kernel");
  params.addCoupledVar("antiferrodis_A_z", "variable antiferrodis_A_z coupled into this kernel");
  params.addCoupledVar("e_xx", "variable e_xx coupled into this kernel");
  params.addCoupledVar("e_yy", "variable e_yy coupled into this kernel");
  params.addCoupledVar("e_zz", "variable e_zz coupled into this kernel");
  params.addCoupledVar("e_xy", "variable e_xy coupled into this kernel");
  params.addCoupledVar("e_yz", "variable e_yz coupled into this kernel");
  params.addCoupledVar("e_zx", "variable e_zx coupled into this kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2.0 for z)");
  params.addParam<Real>("G", 1.0, "the time constant");
  params.addParam<Real>("r11", 1.0, "the first coefficient");
  params.addParam<Real>("r12", 1.0, "the second coefficient");
  params.addParam<Real>("r44", 1.0, "the third coefficient");
  return params;
}

ScalarRotostrictiveEnergy::ScalarRotostrictiveEnergy(const InputParameters & parameters)
  : ODEKernel(parameters),
   _component(getParam<unsigned int>("component")),
    _antiferrodis_A_x_var(coupledScalar("antiferrodis_A_x")),
    _antiferrodis_A_y_var(coupledScalar("antiferrodis_A_x")),
    _antiferrodis_A_z_var(coupledScalar("antiferrodis_A_x")),
    _e_xx_var(coupledScalar("e_xx")),
    _e_yy_var(coupledScalar("e_yy")),
    _e_zz_var(coupledScalar("e_zz")),
    _e_xy_var(coupledScalar("e_xy")),
    _e_yz_var(coupledScalar("e_yz")),
    _e_zx_var(coupledScalar("e_zx")),
    _antiferrodis_A_x(coupledScalarValue("antiferrodis_A_x")),
    _antiferrodis_A_y(coupledScalarValue("antiferrodis_A_y")),
    _antiferrodis_A_z(coupledScalarValue("antiferrodis_A_z")),
    _e_xx(coupledScalarValue("e_xx")),
    _e_yy(coupledScalarValue("e_yy")),
    _e_zz(coupledScalarValue("e_zz")),
    _e_xy(coupledScalarValue("e_xy")),
    _e_yz(coupledScalarValue("e_yz")),
    _e_zx(coupledScalarValue("e_zx")),
    _G(getParam<Real>("G")),
    _r11(getParam<Real>("r11")),
    _r12(getParam<Real>("r12")),
    _r44(getParam<Real>("r44"))
{
}

Real
ScalarRotostrictiveEnergy::computeQpResidual()
{
  if (_component == 0)
  {
    return _G*(2*_e_xx[_i]*_antiferrodis_A_x[_i]*_r11 + (2*_e_yy[_i]*_antiferrodis_A_x[_i] + 2*_e_zz[_i]*_antiferrodis_A_x[_i])*_r12 + (2*_e_zx[_i]*_antiferrodis_A_y[_i] + 2*_e_yz[_i]*_antiferrodis_A_z[_i])*_r44);
  }
  else if (_component == 1)
  {
    return _G*(2*_e_yy[_i]*_antiferrodis_A_y[_i]*_r11 + (2*_e_xx[_i]*_antiferrodis_A_y[_i] + 2*_e_zz[_i]*_antiferrodis_A_y[_i])*_r12 + (2*_e_zx[_i]*_antiferrodis_A_x[_i] + 2*_e_xy[_i]*_antiferrodis_A_z[_i])*_r44);
  }
  else if (_component == 2)
  {
    return _G*(2*_e_zz[_i]*_antiferrodis_A_z[_i]*_r11 + (2*_e_xx[_i]*_antiferrodis_A_z[_i] + 2*_e_yy[_i]*_antiferrodis_A_z[_i])*_r12 + (2*_e_yz[_i]*_antiferrodis_A_x[_i] + 2*_e_xy[_i]*_antiferrodis_A_y[_i])*_r44);
  }
  else if (_component == 3)
  {
    return _G*(Utility::pow<2>(_antiferrodis_A_x[_i])*_r11 + (Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_r12);
  }
  else if (_component == 4)
  {
    return _G*(Utility::pow<2>(_antiferrodis_A_y[_i])*_r11 + (Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_r12);
  }
  else if (_component == 5)
  {
    return _G*(Utility::pow<2>(_antiferrodis_A_z[_i])*_r11 + (Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]))*_r12);
  }
  else if (_component == 6)
  {
    return 2*_G*_antiferrodis_A_y[_i]*_antiferrodis_A_z[_i]*_r44;
  }
  else if (_component == 7)
  {
    return 2*_G*_antiferrodis_A_x[_i]*_antiferrodis_A_z[_i]*_r44;
  }
  else if (_component == 8)
  {
    return 2*_G*_antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]*_r44;
  }
  else
    return 0.0;
}

Real
ScalarRotostrictiveEnergy::computeQpJacobian()
{
  if (_component == 0)
  {
    return _G*(2*_e_xx[_i]*_r11 + (2*_e_yy[_i] + 2*_e_zz[_i])*_r12);
  }
  else if (_component == 1)
  {
    return _G*(2*_e_yy[_i]*_r11 + (2*_e_xx[_i] + 2*_e_zz[_i])*_r12);
  }
  else if (_component == 2)
  {
    return _G*(2*_e_zz[_i]*_r11 + (2*_e_xx[_i] + 2*_e_yy[_i])*_r12);
  }
  else if (_component == 3)
  {
    return 0.0;
  }
  else if (_component == 4)
  {
    return 0.0;
  }
  else if (_component == 5)
  {
    return 0.0;
  }
  else if (_component == 6)
  {
    return 0.0;
  }
  else if (_component == 7)
  {
    return 0.0;
  }
  else if (_component == 8)
  {
    return 0.0;
  }
  else
    return 0.0;
}

Real
ScalarRotostrictiveEnergy::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  if (_component == 0)
  {
   if (jvar == _antiferrodis_A_y_var)
   {
     return 2*_e_zx[_i]*_G*_r44;
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return 2*_e_yz[_i]*_G*_r44;
   }
   else if (jvar == _e_xx_var)
   {
     return 2*_G*_antiferrodis_A_x[_i]*_r11;
   }
   else if (jvar == _e_yy_var)
   {
     return 2*_G*_antiferrodis_A_x[_i]*_r12;
   }
   else if (jvar == _e_zz_var)
   {
     return 2*_G*_antiferrodis_A_x[_i]*_r12;
   }
   else
   {
     return 0.0;
   }
  }
  else if (_component == 1)
  {
   if (jvar == _antiferrodis_A_x_var)
   {
     return 2*_e_zx[_i]*_G*_r44;
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return 2*_e_xy[_i]*_G*_r44;
   }
   else if (jvar == _e_xx_var)
   {
     return 2*_G*_antiferrodis_A_y[_i]*_r12;
   }
   else if (jvar == _e_yy_var)
   {
     return 2*_G*_antiferrodis_A_y[_i]*_r11;
   }
   else if (jvar == _e_zz_var)
   {
     return 2*_G*_antiferrodis_A_y[_i]*_r12;
   }
   else
   {
     return 0.0;
   }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return 2*_e_yz[_i]*_G*_r44;
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return 2*_e_xy[_i]*_G*_r44;
    }
    else if (jvar == _e_xx_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r12;
    }
    else if (jvar == _e_yy_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r12;
    }
    else if (jvar == _e_zz_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r11;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 3)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return 2*_G*_antiferrodis_A_x[_i]*_r11;
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return 2*_G*_antiferrodis_A_y[_i]*_r12;
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r12;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 4)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return 2*_G*_antiferrodis_A_x[_i]*_r12;
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return 2*_G*_antiferrodis_A_y[_i]*_r11;
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r12;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 5)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r12;
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r12;
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r11;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 6)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return 0.0;
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r44;
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return 2*_G*_antiferrodis_A_y[_i]*_r44;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 7)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return 2*_G*_antiferrodis_A_z[_i]*_r44;
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return 0.0;
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return 2*_G*_antiferrodis_A_x[_i]*_r44;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 8)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return 2*_G*_antiferrodis_A_y[_i]*_r44;
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return 2*_G*_antiferrodis_A_x[_i]*_r44;
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return 0.0;
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
