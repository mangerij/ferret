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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ScalarRotopolarEnergy.h"

registerMooseObject("FerretApp", ScalarRotopolarEnergy);

InputParameters
ScalarRotopolarEnergy::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addCoupledVar("polar_x", "variable polar_x coupled into this kernel");
  params.addCoupledVar("polar_y", "variable polar_y coupled into this kernel");
  params.addCoupledVar("polar_z", "variable polar_z coupled into this kernel");
  params.addCoupledVar("antiferrodis_A_x", "variable antiferrodis_A_x coupled into this kernel");
  params.addCoupledVar("antiferrodis_A_y", "variable antiferrodis_A_y coupled into this kernel");
  params.addCoupledVar("antiferrodis_A_z", "variable antiferrodis_A_z coupled into this kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for _polar_x[_i], 1 for _polar_y[_i], 2 for _polar_z[_i], 3 for _antiferrodis_A_x[_i], 4 for _antiferrodis_A_y[_i], 5 for _antiferrodis_A_z[_i])");
  params.addParam<Real>("G", 1.0, "the time constant");
  params.addParam<Real>("Bpr", 1.0, "the first coefficient");
  params.addParam<Real>("Cpr", 1.0, "the second coefficient");
  params.addParam<Real>("Cprp", 1.0, "the third coefficient");
  return params;
}

ScalarRotopolarEnergy::ScalarRotopolarEnergy(const InputParameters & parameters)
  : ODEKernel(parameters),
   _component(getParam<unsigned int>("component")),
    _polar_x_var(coupledScalar("polar_x")),
    _polar_y_var(coupledScalar("polar_x")),
    _polar_z_var(coupledScalar("polar_x")),
    _antiferrodis_A_x_var(coupledScalar("antiferrodis_A_x")),
    _antiferrodis_A_y_var(coupledScalar("antiferrodis_A_y")),
    _antiferrodis_A_z_var(coupledScalar("antiferrodis_A_z")),
    _polar_x(coupledScalarValue("polar_x")),
    _polar_y(coupledScalarValue("polar_y")),
    _polar_z(coupledScalarValue("polar_z")),
    _antiferrodis_A_x(coupledScalarValue("antiferrodis_A_x")),
    _antiferrodis_A_y(coupledScalarValue("antiferrodis_A_y")),
    _antiferrodis_A_z(coupledScalarValue("antiferrodis_A_z")),
    _G(getParam<Real>("G")),
    _Bpr(getParam<Real>("Bpr")),
    _Cpr(getParam<Real>("Cpr")),
    _Cprp(getParam<Real>("Cprp"))
{
}

Real
ScalarRotopolarEnergy::computeQpResidual()
{
  if (_component == 0)
  {
    return _G*(2*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_Bpr*_polar_x[_i] + 2*Utility::pow<2>(_antiferrodis_A_x[_i])*_Cpr*_polar_x[_i] + _Cprp*(_antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]*_polar_y[_i] + _antiferrodis_A_x[_i]*_antiferrodis_A_z[_i]*_polar_z[_i]));
  }
  else if (_component == 1)
  {
    return _G*(2*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_Bpr*_polar_y[_i] + 2*Utility::pow<2>(_antiferrodis_A_y[_i])*_Cpr*_polar_y[_i] + _Cprp*(_antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]*_polar_x[_i] + _antiferrodis_A_y[_i]*_antiferrodis_A_z[_i]*_polar_z[_i]));
  }
  else if (_component == 2)
  {
    return _G*(_Cprp*(_antiferrodis_A_x[_i]*_antiferrodis_A_z[_i]*_polar_x[_i] + _antiferrodis_A_y[_i]*_antiferrodis_A_z[_i]*_polar_y[_i]) + 2*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_Bpr*_polar_z[_i] + 2*Utility::pow<2>(_antiferrodis_A_z[_i])*_Cpr*_polar_z[_i]);
  }
  else if (_component == 3)
  {
    return _G*(2*_antiferrodis_A_x[_i]*_Cpr*Utility::pow<2>(_polar_x[_i]) + _Cprp*(_antiferrodis_A_y[_i]*_polar_x[_i]*_polar_y[_i] + _antiferrodis_A_z[_i]*_polar_x[_i]*_polar_z[_i]) + 2*_antiferrodis_A_x[_i]*_Bpr*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 4)
  {
    return _G*(2*_antiferrodis_A_y[_i]*_Cpr*Utility::pow<2>(_polar_y[_i]) + _Cprp*(_antiferrodis_A_x[_i]*_polar_x[_i]*_polar_y[_i] + _antiferrodis_A_z[_i]*_polar_y[_i]*_polar_z[_i]) + 2*_antiferrodis_A_y[_i]*_Bpr*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 5)
  {
    return _G*(2*_antiferrodis_A_z[_i]*_Cpr*Utility::pow<2>(_polar_z[_i]) + _Cprp*(_antiferrodis_A_x[_i]*_polar_x[_i]*_polar_z[_i] + _antiferrodis_A_y[_i]*_polar_y[_i]*_polar_z[_i]) + 2*_antiferrodis_A_z[_i]*_Bpr*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else
    return 0.0;
}

Real
ScalarRotopolarEnergy::computeQpJacobian()
{
  if (_component == 0)
  {
    return (2*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_Bpr + 2*Utility::pow<2>(_antiferrodis_A_x[_i])*_Cpr)*_G;
  }
  else if (_component == 1)
  {
    return (2*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_Bpr + 2*Utility::pow<2>(_antiferrodis_A_y[_i])*_Cpr)*_G;
  }
  else if (_component == 2)
  {
    return (2*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i]))*_Bpr + 2*Utility::pow<2>(_antiferrodis_A_z[_i])*_Cpr)*_G;
  }
  else if (_component == 3)
  {
    return _G*(2*_Cpr*Utility::pow<2>(_polar_x[_i]) + 2*_Bpr*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 4)
  {
    return _G*(2*_Cpr*Utility::pow<2>(_polar_x[_i]) + 2*_Bpr*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 5)
  {
    return _G*(2*_Cpr*Utility::pow<2>(_polar_x[_i]) + 2*_Bpr*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else
    return 0.0;
}

Real
ScalarRotopolarEnergy::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  if (_component == 0)
  {
   if (jvar == _polar_y_var)
   {
     return _antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]*_Cprp*_G;
   }
   else if (jvar == _polar_z_var)
   {
     return _antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]*_Cprp*_G;
   }
   else if (jvar == _antiferrodis_A_x_var)
   {
     return _G*(4*_antiferrodis_A_x[_i]*_Bpr*_polar_x[_i] + 4*_antiferrodis_A_x[_i]*_Cpr*_polar_x[_i] + _Cprp*(_antiferrodis_A_y[_i]*_polar_y[_i] + _antiferrodis_A_z[_i]*_polar_z[_i]));
   }
   else if (jvar == _antiferrodis_A_y_var)
   {
     return _G*(4*_antiferrodis_A_y[_i]*_Bpr*_polar_x[_i] + _antiferrodis_A_x[_i]*_Cprp*_polar_y[_i]);
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return _G*(4*_antiferrodis_A_z[_i]*_Bpr*_polar_x[_i] + _antiferrodis_A_x[_i]*_Cprp*_polar_z[_i]);
   }
   else
     return 0.0;
  }
  else if (_component == 1)
  {
   if (jvar == _polar_x_var)
   {
     return _antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]*_Cprp*_G;
   }
   else if (jvar == _polar_z_var)
   {
     return _antiferrodis_A_y[_i]*_antiferrodis_A_z[_i]*_Cprp*_G;
   }
   else if (jvar == _antiferrodis_A_x_var)
   {
     return _G*(_antiferrodis_A_y[_i]*_Cprp*_polar_x[_i] + 4*_antiferrodis_A_x[_i]*_Bpr*_polar_y[_i]);
   }
   else if (jvar == _antiferrodis_A_y_var)
   {
     return _G*(4*_antiferrodis_A_y[_i]*_Bpr*_polar_y[_i] + 4*_antiferrodis_A_y[_i]*_Cpr*_polar_y[_i] + _Cprp*(_antiferrodis_A_x[_i]*_polar_x[_i] + _antiferrodis_A_z[_i]*_polar_z[_i]));
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return _G*(4*_antiferrodis_A_z[_i]*_Bpr*_polar_y[_i] + _antiferrodis_A_y[_i]*_Cprp*_polar_z[_i]);
   }
   else
     return 0.0;
  }
  else if (_component == 2)
  {
   if (jvar == _polar_x_var)
   {
     return _antiferrodis_A_x[_i]*_antiferrodis_A_z[_i]*_Cprp*_G;
   }
   else if (jvar == _polar_y_var)
   {
     return _antiferrodis_A_y[_i]*_antiferrodis_A_z[_i]*_Cprp*_G;
   }
   else if (jvar == _antiferrodis_A_x_var)
   {
     return _G*(_antiferrodis_A_z[_i]*_Cprp*_polar_x[_i] + 4*_antiferrodis_A_x[_i]*_Bpr*_polar_z[_i]);
   }
   else if (jvar == _antiferrodis_A_y_var)
   {
     return _G*(_antiferrodis_A_z[_i]*_Cprp*_polar_y[_i] + 4*_antiferrodis_A_y[_i]*_Bpr*_polar_z[_i]);
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return _G*(_Cprp*(_antiferrodis_A_x[_i]*_polar_x[_i] + _antiferrodis_A_y[_i]*_polar_y[_i]) + 4*_antiferrodis_A_z[_i]*_Bpr*_polar_z[_i] + 4*_antiferrodis_A_z[_i]*_Cpr*_polar_z[_i]);
   }
   else
     return 0.0;
  }
  else if (_component == 3)
  {
   if (jvar == _antiferrodis_A_y_var)
   {
     return _Cprp*_G*_polar_x[_i]*_polar_y[_i];
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return _Cprp*_G*_polar_x[_i]*_polar_z[_i];
   }
   else if (jvar == _polar_x_var)
   {
     return _G*(4*_antiferrodis_A_x[_i]*_Bpr*_polar_x[_i] + 4*_antiferrodis_A_x[_i]*_Cpr*_polar_x[_i] + _Cprp*(_antiferrodis_A_y[_i]*_polar_y[_i] + _antiferrodis_A_z[_i]*_polar_z[_i]));
   }
   else if (jvar == _polar_y_var)
   {
     return _G*(_antiferrodis_A_y[_i]*_Cprp*_polar_x[_i] + 4*_antiferrodis_A_x[_i]*_Bpr*_polar_y[_i]);
   }
   else if (jvar == _polar_z_var)
   {
     return _G*(_antiferrodis_A_z[_i]*_Cprp*_polar_x[_i] + 4*_antiferrodis_A_x[_i]*_Bpr*_polar_z[_i]);
   }
   else
     return 0.0;
  }
  else if (_component == 4)
  {
   if (jvar == _antiferrodis_A_x_var)
   {
     return _Cprp*_G*_polar_x[_i]*_polar_y[_i];
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return _Cprp*_G*_polar_y[_i]*_polar_z[_i];
   }
   else if (jvar == _polar_x_var)
   {
     return _G*(4*_antiferrodis_A_y[_i]*_Bpr*_polar_x[_i] + _antiferrodis_A_x[_i]*_Cprp*_polar_y[_i]);
   }
   else if (jvar == _polar_y_var)
   {
     return _G*(4*_antiferrodis_A_y[_i]*_Bpr*_polar_y[_i] + 4*_antiferrodis_A_y[_i]*_Cpr*_polar_y[_i] + _Cprp*(_antiferrodis_A_x[_i]*_polar_x[_i] + _antiferrodis_A_z[_i]*_polar_z[_i]));
   }
   else if (jvar == _polar_z_var)
   {
     return _G*(_antiferrodis_A_z[_i]*_Cprp*_polar_y[_i] + 4*_antiferrodis_A_y[_i]*_Bpr*_polar_z[_i]);
   }
   else
     return 0.0;
  }
  else if (_component == 5)
  {
   if (jvar == _antiferrodis_A_x_var)
   {
     return _Cprp*_G*_polar_x[_i]*_polar_z[_i];
   }
   else if (jvar == _antiferrodis_A_y_var)
   {
     return _Cprp*_G*_polar_y[_i]*_polar_z[_i];
   }
   else if (jvar == _polar_x_var)
   {
     return _G*(4*_antiferrodis_A_z[_i]*_Bpr*_polar_x[_i] + _antiferrodis_A_x[_i]*_Cprp*_polar_z[_i]);
   }
   else if (jvar == _polar_y_var)
   {
     return _G*(4*_antiferrodis_A_z[_i]*_Bpr*_polar_y[_i] + _antiferrodis_A_y[_i]*_Cprp*_polar_z[_i]);
   }
   else if (jvar == _polar_z_var)
   {
     return _G*(_Cprp*(_antiferrodis_A_x[_i]*_polar_x[_i] + _antiferrodis_A_y[_i]*_polar_y[_i]) + 4*_antiferrodis_A_z[_i]*_Bpr*_polar_z[_i] + 4*_antiferrodis_A_z[_i]*_Cpr*_polar_z[_i]);
   }
   else
     return 0.0;
  }
  else
    return 0.0;
}
