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

#include "ScalarElasticEnergy.h"

registerMooseObject("FerretApp", ScalarElasticEnergy);

template <>
InputParameters
validParams<ScalarElasticEnergy>()
{
  InputParameters params = validParams<ODEKernel>();
  params.addCoupledVar("e_xx", "variable e_xx coupled into this kernel");
  params.addCoupledVar("e_yy", "variable e_yy coupled into this kernel");
  params.addCoupledVar("e_zz", "variable e_zz coupled into this kernel");
  params.addCoupledVar("e_xy", "variable e_xy coupled into this kernel");
  params.addCoupledVar("e_yz", "variable e_yz coupled into this kernel");
  params.addCoupledVar("e_zx", "variable e_zx coupled into this kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for xx, 1 for yy, 2 for zz, 3 for xy, 4 for yz, 5 for zx)");
  params.addParam<Real>("G", 1.0, "the time constant");
  params.addParam<Real>("C11", 1.0, "the first coefficient");
  params.addParam<Real>("C12", 1.0, "the second coefficient");
  params.addParam<Real>("C44", 1.0, "the third coefficient");
  return params;
}

ScalarElasticEnergy::ScalarElasticEnergy(const InputParameters & parameters)
  : ODEKernel(parameters),
   _component(getParam<unsigned int>("component")),
    _e_xx_var(coupledScalar("e_xx")),
    _e_yy_var(coupledScalar("e_yy")),
    _e_zz_var(coupledScalar("e_zz")),
    _e_xy_var(coupledScalar("e_xy")),
    _e_yz_var(coupledScalar("e_yz")),
    _e_zx_var(coupledScalar("e_zx")),
    _e_xx(coupledScalarValue("e_xx")),
    _e_yy(coupledScalarValue("e_yy")),
    _e_zz(coupledScalarValue("e_zz")),
    _e_xy(coupledScalarValue("e_xy")),
    _e_yz(coupledScalarValue("e_yz")),
    _e_zx(coupledScalarValue("e_zx")),
    _G(getParam<Real>("G")),
    _C11(getParam<Real>("C11")),
    _C12(getParam<Real>("C12")),
    _C44(getParam<Real>("C44"))
{
}

Real
ScalarElasticEnergy::computeQpResidual()
{
  if (_component == 0)
  {
    return (_C11*_e_xx[_i] + _C12*(_e_yy[_i] + _e_zz[_i]))*_G;
  }
  else if (_component == 1)
  {
    return (_C11*_e_yy[_i] + _C12*(_e_xx[_i] + _e_zz[_i]))*_G;
  }
  else if (_component == 2)
  {
    return (_C12*(_e_xx[_i] + _e_yy[_i]) + _C11*_e_zz[_i])*_G;
  }
  else if (_component == 3)
  {
    return 4*_C44*_e_xy[_i]*_G;
  }
  else if (_component == 4)
  {
    return 4*_C44*_e_yz[_i]*_G;
  }
  else if (_component == 5)
  {
    return 4*_C44*_e_zx[_i]*_G;
  }
  else
    return 0.0;
}

Real
ScalarElasticEnergy::computeQpJacobian()
{
  if (_component == 0)
  {
    return _C11*_G;
  }
  else if (_component == 1)
  {
    return _C11*_G;
  }
  else if (_component == 2)
  {
    return _C11*_G;
  }
  else if (_component == 3)
  {
    return 4*_C44*_G;
  }
  else if (_component == 4)
  {
    return 4*_C44*_G;
  }
  else if (_component == 5)
  {
    return 4*_C44*_G;
  }
  else
    return 0.0;
}

Real
ScalarElasticEnergy::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  if (_component == 0)
  {
   if (jvar == _e_yy_var)
   {
     return _C12*_G;
   }
   else if (jvar == _e_zz_var)
   {
     return _C12*_G;
   }
   else
   {
     return 0.0;
   }
  }
  else if (_component == 1)
  {
   if (jvar == _e_xx_var)
   {
     return _C12*_G;
   }
   else if (jvar == _e_zz_var)
   {
     return _C12*_G;
   }
   else
   {
     return 0.0;
   }
  }
  else if (_component == 2)
  {
   if (jvar == _e_xx_var)
   {
     return _C12*_G;
   }
   else if (jvar == _e_yy_var)
   {
     return _C12*_G;
   }
   else
   {
     return 0.0;
   }
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
  else
    return 0.0;
}
