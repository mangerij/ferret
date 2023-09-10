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

   You should have received a co_polar_y[_i] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ScalarElectrostrictiveEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ScalarElectrostrictiveEnergy);

InputParameters
ScalarElectrostrictiveEnergy::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addCoupledVar("polar_x", "variable polar_x coupled into this kernel");
  params.addCoupledVar("polar_y", "variable polar_y coupled into this kernel");
  params.addCoupledVar("polar_z", "variable polar_z coupled into this kernel");
  params.addCoupledVar("e_xx", "variable e_xx coupled into this kernel");
  params.addCoupledVar("e_yy", "variable e_yy coupled into this kernel");
  params.addCoupledVar("e_zz", "variable e_zz coupled into this kernel");
  params.addCoupledVar("e_xy", "variable e_xy coupled into this kernel");
  params.addCoupledVar("e_yz", "variable e_yz coupled into this kernel");
  params.addCoupledVar("e_zx", "variable e_zx coupled into this kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the variable this kernel acts in. (0 for _polar_x[_i], 1 for _polar_y[_i], 2 for _polar_z[_i], 3 _e_xx[_i], 4 _e_yy[_i],5 _e_zz[_i],6 _e_xy[_i],7 _e_yz[_i],8 ezy)");
  params.addParam<Real>("G", 1.0, "the time constant");
  params.addParam<Real>("q11", 1.0, "the first coefficient");
  params.addParam<Real>("q12", 1.0, "the second coefficient");
  params.addParam<Real>("q44", 1.0, "the third coefficient");
  return params;
}

ScalarElectrostrictiveEnergy::ScalarElectrostrictiveEnergy(const InputParameters & parameters)
  : ODEKernel(parameters),
   _component(getParam<unsigned int>("component")),
    _polar_x_var(coupledScalar("polar_x")),
    _polar_y_var(coupledScalar("polar_x")),
    _polar_z_var(coupledScalar("polar_x")),
    _e_xx_var(coupledScalar("e_xx")),
    _e_yy_var(coupledScalar("e_yy")),
    _e_zz_var(coupledScalar("e_zz")),
    _e_xy_var(coupledScalar("e_xy")),
    _e_yz_var(coupledScalar("e_yz")),
    _e_zx_var(coupledScalar("e_zx")),
    _polar_x(coupledScalarValue("polar_x")),
    _polar_y(coupledScalarValue("polar_y")),
    _polar_z(coupledScalarValue("polar_z")),
    _e_xx(coupledScalarValue("e_xx")),
    _e_yy(coupledScalarValue("e_yy")),
    _e_zz(coupledScalarValue("e_zz")),
    _e_xy(coupledScalarValue("e_xy")),
    _e_yz(coupledScalarValue("e_yz")),
    _e_zx(coupledScalarValue("e_zx")),
    _G(getParam<Real>("G")),
    _q11(getParam<Real>("q11")),
    _q12(getParam<Real>("q12")),
    _q44(getParam<Real>("q44"))
{
}

Real
ScalarElectrostrictiveEnergy::computeQpResidual()
{
  if (_component == 0)
  {
    return _G*(2*_e_xx[_i]*_polar_x[_i]*_q11 + (2*_e_yy[_i]*_polar_x[_i] + 2*_e_zz[_i]*_polar_x[_i])*_q12 + (2*_e_zx[_i]*_polar_y[_i] + 2*_e_yz[_i]*_polar_z[_i])*_q44);
  }
  else if (_component == 1)
  {
    return _G*(2*_e_yy[_i]*_polar_y[_i]*_q11 + (2*_e_xx[_i]*_polar_y[_i] + 2*_e_zz[_i]*_polar_y[_i])*_q12 + (2*_e_zx[_i]*_polar_x[_i] + 2*_e_xy[_i]*_polar_z[_i])*_q44);
  }
  else if (_component == 2)
  {
    return _G*(2*_e_zz[_i]*_polar_z[_i]*_q11 + (2*_e_xx[_i]*_polar_z[_i] + 2*_e_yy[_i]*_polar_z[_i])*_q12 + (2*_e_yz[_i]*_polar_x[_i] + 2*_e_xy[_i]*_polar_y[_i])*_q44);
  }
  else if (_component == 3)
  {
    return _G*(Utility::pow<2>(_polar_x[_i])*_q11 + (Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i]))*_q12);
  }
  else if (_component == 4)
  {
    return _G*(Utility::pow<2>(_polar_y[_i])*_q11 + (Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_z[_i]))*_q12);
  }
  else if (_component == 5)
  {
    return _G*(Utility::pow<2>(_polar_z[_i])*_q11 + (Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]))*_q12);
  }
  else if (_component == 6)
  {
    return 2*_G*_polar_y[_i]*_polar_z[_i]*_q44;
  }
  else if (_component == 7)
  {
    return 2*_G*_polar_x[_i]*_polar_z[_i]*_q44;
  }
  else if (_component == 8)
  {
    return 2*_G*_polar_x[_i]*_polar_y[_i]*_q44;
  }
  else
    return 0.0;
}

Real
ScalarElectrostrictiveEnergy::computeQpJacobian()
{
  if (_component == 0)
  {
    return _G*(2*_e_xx[_i]*_q11 + (2*_e_yy[_i] + 2*_e_zz[_i])*_q12);
  }
  else if (_component == 1)
  {
    return _G*(2*_e_yy[_i]*_q11 + (2*_e_xx[_i] + 2*_e_zz[_i])*_q12);
  }
  else if (_component == 2)
  {
    return _G*(2*_e_zz[_i]*_q11 + (2*_e_xx[_i] + 2*_e_yy[_i])*_q12);
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
ScalarElectrostrictiveEnergy::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  if (_component == 0)
  {
   if (jvar == _polar_y_var)
   {
     return 2*_e_zx[_i]*_G*_q44;
   }
   else if (jvar == _polar_z_var)
   {
     return 2*_e_yz[_i]*_G*_q44;
   }
   else if (jvar == _e_xx_var)
   {
     return 2*_G*_polar_x[_i]*_q11;
   }
   else if (jvar == _e_yy_var)
   {
     return 2*_G*_polar_x[_i]*_q12;
   }
   else if (jvar == _e_zz_var)
   {
     return 2*_G*_polar_x[_i]*_q12;
   }
   else
   {
     return 0.0;
   }
  }
  else if (_component == 1)
  {
   if (jvar == _polar_x_var)
   {
     return 2*_e_zx[_i]*_G*_q44;
   }
   else if (jvar == _polar_z_var)
   {
     return 2*_e_xy[_i]*_G*_q44;
   }
   else if (jvar == _e_xx_var)
   {
     return 2*_G*_polar_y[_i]*_q12;
   }
   else if (jvar == _e_yy_var)
   {
     return 2*_G*_polar_y[_i]*_q11;
   }
   else if (jvar == _e_zz_var)
   {
     return 2*_G*_polar_y[_i]*_q12;
   }
   else
   {
     return 0.0;
   }
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return 2*_e_yz[_i]*_G*_q44;
    }
    else if (jvar == _polar_y_var)
    {
      return 2*_e_xy[_i]*_G*_q44;
    }
    else if (jvar == _e_xx_var)
    {
      return 2*_G*_polar_z[_i]*_q12;
    }
    else if (jvar == _e_yy_var)
    {
      return 2*_G*_polar_z[_i]*_q12;
    }
    else if (jvar == _e_zz_var)
    {
      return 2*_G*_polar_z[_i]*_q11;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 3)
  {
    if (jvar == _polar_x_var)
    {
      return 2*_G*_polar_x[_i]*_q11;
    }
    else if (jvar == _polar_y_var)
    {
      return 2*_G*_polar_y[_i]*_q12;
    }
    else if (jvar == _polar_z_var)
    {
      return 2*_G*_polar_z[_i]*_q12;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 4)
  {
    if (jvar == _polar_x_var)
    {
      return 2*_G*_polar_x[_i]*_q12;
    }
    else if (jvar == _polar_y_var)
    {
      return 2*_G*_polar_y[_i]*_q11;
    }
    else if (jvar == _polar_z_var)
    {
      return 2*_G*_polar_z[_i]*_q12;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 5)
  {
    if (jvar == _polar_x_var)
    {
      return 2*_G*_polar_z[_i]*_q12;
    }
    else if (jvar == _polar_y_var)
    {
      return 2*_G*_polar_z[_i]*_q12;
    }
    else if (jvar == _polar_z_var)
    {
      return 2*_G*_polar_z[_i]*_q11;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 6)
  {
    if (jvar == _polar_x_var)
    {
      return 0.0;
    }
    else if (jvar == _polar_y_var)
    {
      return 2*_G*_polar_z[_i]*_q44;
    }
    else if (jvar == _polar_z_var)
    {
      return 2*_G*_polar_y[_i]*_q44;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 7)
  {
    if (jvar == _polar_x_var)
    {
      return 2*_G*_polar_z[_i]*_q44;
    }
    else if (jvar == _polar_y_var)
    {
      return 0.0;
    }
    else if (jvar == _polar_z_var)
    {
      return 2*_G*_polar_x[_i]*_q44;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 8)
  {
    if (jvar == _polar_x_var)
    {
      return 2*_G*_polar_y[_i]*_q44;
    }
    else if (jvar == _polar_y_var)
    {
      return 2*_G*_polar_x[_i]*_q44;
    }
    else if (jvar == _polar_z_var)
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
