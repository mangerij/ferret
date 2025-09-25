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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ScalarBulkEnergyP.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ScalarBulkEnergyP);

InputParameters
ScalarBulkEnergyP::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addCoupledVar("polar_x", "variable polar_x coupled into this kernel");
  params.addCoupledVar("polar_y", "variable polar_y coupled into this kernel");
  params.addCoupledVar("polar_z", "variable polar_z coupled into this kernel");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2.0 for z)");
  params.addParam<Real>("G", 1.0, "the time constant");
  params.addParam<Real>("Ap", 1.0, "the first coefficient");
  params.addParam<Real>("Bp", 1.0, "the second coefficient");
  params.addParam<Real>("Cp", 1.0, "the third coefficient");
  return params;
}

ScalarBulkEnergyP::ScalarBulkEnergyP(const InputParameters & parameters)
  : ODEKernel(parameters),
   _component(getParam<unsigned int>("component")),
    _polar_x_var(coupledScalar("polar_x")),
    _polar_y_var(coupledScalar("polar_x")),
    _polar_z_var(coupledScalar("polar_x")),
    _polar_x(coupledScalarValue("polar_x")),
    _polar_y(coupledScalarValue("polar_y")),
    _polar_z(coupledScalarValue("polar_z")),
    _G(getParam<Real>("G")),
    _Ap(getParam<Real>("Ap")),
    _Bp(getParam<Real>("Bp")),
    _Cp(getParam<Real>("Cp"))
{
}

Real
ScalarBulkEnergyP::computeQpResidual()
{
  if (_component == 0)
  {
    return _G*(2*_Ap*_polar_x[_i] + 4*_Bp*_polar_x[_i]*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])) + _Cp*(2*_polar_x[_i]*Utility::pow<2>(_polar_y[_i]) + 2*_polar_x[_i]*Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 1)
  {
    return _G*(2*_Ap*_polar_y[_i] + 4*_Bp*_polar_y[_i]*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])) + _Cp*(2*Utility::pow<2>(_polar_x[_i])*_polar_y[_i] + 2*_polar_y[_i]*Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 2)
  {
    return _G*(2*_Ap*_polar_z[_i] + _Cp*(2*Utility::pow<2>(_polar_x[_i])*_polar_z[_i] + 2*Utility::pow<2>(_polar_y[_i])*_polar_z[_i]) + 4*_Bp*_polar_z[_i]*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else
    return 0.0;
}

Real
ScalarBulkEnergyP::computeQpJacobian()
{
  if (_component == 0)
  {
    return _G*(2*_Ap + 8*_Bp*Utility::pow<2>(_polar_x[_i]) + 4*_Bp*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])) + _Cp*(2*Utility::pow<2>(_polar_y[_i]) + 2*Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 1)
  {
    return _G*(2*_Ap + 8*_Bp*Utility::pow<2>(_polar_y[_i]) + 4*_Bp*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])) + _Cp*(2*Utility::pow<2>(_polar_x[_i]) + 2*Utility::pow<2>(_polar_z[_i])));
  }
  else if (_component == 2)
  {
    return _G*(2*_Ap + _Cp*(2*Utility::pow<2>(_polar_x[_i]) + 2*Utility::pow<2>(_polar_y[_i])) + 8*_Bp*Utility::pow<2>(_polar_z[_i]) + 4*_Bp*(Utility::pow<2>(_polar_x[_i]) + Utility::pow<2>(_polar_y[_i]) + Utility::pow<2>(_polar_z[_i])));
  }
  else
    return 0.0;
}

Real
ScalarBulkEnergyP::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _G*(8*_Bp*_polar_x[_i]*_polar_y[_i] + 4*_Cp*_polar_x[_i]*_polar_y[_i]);
    }
    else if (jvar == _polar_z_var)
    {
      return _G*(8*_Bp*_polar_x[_i]*_polar_z[_i] + 4*_Cp*_polar_x[_i]*_polar_z[_i]);
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
      return _G*(8*_Bp*_polar_x[_i]*_polar_y[_i] + 4*_Cp*_polar_x[_i]*_polar_y[_i]);
    }
    else if (jvar == _polar_z_var)
    {
      return _G*(8*_Bp*_polar_y[_i]*_polar_z[_i] + 4*_Cp*_polar_y[_i]*_polar_z[_i]);
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
      return _G*(8*_Bp*_polar_x[_i]*_polar_z[_i] + 4*_Cp*_polar_x[_i]*_polar_z[_i]);
    }
    else if (jvar == _polar_y_var)
    {
      return _G*(8*_Bp*_polar_y[_i]*_polar_z[_i] + 4*_Cp*_polar_y[_i]*_polar_z[_i]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
