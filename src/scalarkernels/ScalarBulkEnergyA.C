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

#include "ScalarBulkEnergyA.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ScalarBulkEnergyA);

InputParameters
ScalarBulkEnergyA::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addCoupledVar("antiferrodis_A_x", "variable antiferrodis_A_x coupled into this kernel");
  params.addCoupledVar("antiferrodis_A_y", "variable antiferrodis_A_y coupled into this kernel");
  params.addCoupledVar("antiferrodis_A_z", "variable antiferrodis_A_z coupled into this kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2.0 for z)");
  params.addParam<Real>("G", 1.0, "the time constant");
  params.addParam<Real>("Ar", 1.0, "the first coefficient");
  params.addParam<Real>("Br", 1.0, "the second coefficient");
  params.addParam<Real>("Cr", 1.0, "the third coefficient");
  return params;
}

ScalarBulkEnergyA::ScalarBulkEnergyA(const InputParameters & parameters)
  : ODEKernel(parameters),
   _component(getParam<unsigned int>("component")),
    _antiferrodis_A_x_var(coupledScalar("antiferrodis_A_x")),
    _antiferrodis_A_y_var(coupledScalar("antiferrodis_A_x")),
    _antiferrodis_A_z_var(coupledScalar("antiferrodis_A_x")),
    _antiferrodis_A_x(coupledScalarValue("antiferrodis_A_x")),
    _antiferrodis_A_y(coupledScalarValue("antiferrodis_A_y")),
    _antiferrodis_A_z(coupledScalarValue("antiferrodis_A_z")),
    _G(getParam<Real>("G")),
    _Ar(getParam<Real>("Ar")),
    _Br(getParam<Real>("Br")),
    _Cr(getParam<Real>("Cr"))
{
}

Real
ScalarBulkEnergyA::computeQpResidual()
{
  if (_component == 0)
  {
    return _G*(2*_Ar*_antiferrodis_A_x[_i] + 4*_Br*_antiferrodis_A_x[_i]*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i])) + _Cr*(2*_antiferrodis_A_x[_i]*Utility::pow<2>(_antiferrodis_A_y[_i]) + 2*_antiferrodis_A_x[_i]*Utility::pow<2>(_antiferrodis_A_z[_i])));
  }
  else if (_component == 1)
  {
    return _G*(2*_Ar*_antiferrodis_A_y[_i] + 4*_Br*_antiferrodis_A_y[_i]*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i])) + _Cr*(2*Utility::pow<2>(_antiferrodis_A_x[_i])*_antiferrodis_A_y[_i] + 2*_antiferrodis_A_y[_i]*Utility::pow<2>(_antiferrodis_A_z[_i])));
  }
  else if (_component == 2)
  {
    return _G*(2*_Ar*_antiferrodis_A_z[_i] + _Cr*(2*Utility::pow<2>(_antiferrodis_A_x[_i])*_antiferrodis_A_z[_i] + 2*Utility::pow<2>(_antiferrodis_A_y[_i])*_antiferrodis_A_z[_i]) + 4*_Br*_antiferrodis_A_z[_i]*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i])));
  }
  else
    return 0.0;
}

Real
ScalarBulkEnergyA::computeQpJacobian()
{
  if (_component == 0)
  {
    return _G*(2*_Ar + 8*_Br*Utility::pow<2>(_antiferrodis_A_x[_i]) + 4*_Br*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i])) + _Cr*(2*Utility::pow<2>(_antiferrodis_A_y[_i]) + 2*Utility::pow<2>(_antiferrodis_A_z[_i])));
  }
  else if (_component == 1)
  {
    return _G*(2*_Ar + 8*_Br*Utility::pow<2>(_antiferrodis_A_y[_i]) + 4*_Br*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i])) + _Cr*(2*Utility::pow<2>(_antiferrodis_A_x[_i]) + 2*Utility::pow<2>(_antiferrodis_A_z[_i])));
  }
  else if (_component == 2)
  {
    return _G*(2*_Ar + _Cr*(2*Utility::pow<2>(_antiferrodis_A_x[_i]) + 2*Utility::pow<2>(_antiferrodis_A_y[_i])) + 8*_Br*Utility::pow<2>(_antiferrodis_A_z[_i]) + 4*_Br*(Utility::pow<2>(_antiferrodis_A_x[_i]) + Utility::pow<2>(_antiferrodis_A_y[_i]) + Utility::pow<2>(_antiferrodis_A_z[_i])));
  }
  else
    return 0.0;
}

Real
ScalarBulkEnergyA::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  if (_component == 0)
  {
   if (jvar == _antiferrodis_A_y_var)
   {
     return _G*(8*_Br*_antiferrodis_A_x[_i]*_antiferrodis_A_y[_i] + 4*_Cr*_antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]);
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return _G*(8*_Br*_antiferrodis_A_x[_i]*_antiferrodis_A_z[_i] + 4*_Cr*_antiferrodis_A_x[_i]*_antiferrodis_A_z[_i]);
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
     return _G*(8*_Br*_antiferrodis_A_x[_i]*_antiferrodis_A_y[_i] + 4*_Cr*_antiferrodis_A_x[_i]*_antiferrodis_A_y[_i]);
   }
   else if (jvar == _antiferrodis_A_z_var)
   {
     return _G*(8*_Br*_antiferrodis_A_y[_i]*_antiferrodis_A_z[_i] + 4*_Cr*_antiferrodis_A_y[_i]*_antiferrodis_A_z[_i]);
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
     return _G*(8*_Br*_antiferrodis_A_x[_i]*_antiferrodis_A_z[_i] + 4*_Cr*_antiferrodis_A_x[_i]*_antiferrodis_A_z[_i]);
   }
   else if (jvar == _antiferrodis_A_y_var)
   {
     return _G*(8*_Br*_antiferrodis_A_y[_i]*_antiferrodis_A_z[_i] + 4*_Cr*_antiferrodis_A_y[_i]*_antiferrodis_A_z[_i]);
   }
   else
   {
     return 0.0;
   }
  }
  else
    return 0.0;
}
