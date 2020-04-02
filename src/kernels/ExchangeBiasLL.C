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

#include "ExchangeBiasLL.h"
#include "IntegratedBC.h"
#include "MooseMesh.h"

#include "libmesh/utility.h"

registerMooseObject("FerretApp", ExchangeBiasLL);

template<>
InputParameters validParams<ExchangeBiasLL>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

ExchangeBiasLL::ExchangeBiasLL(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _g0(getMaterialProperty<Real>("g0")),
  _eEB(getMaterialProperty<Real>("eEB")),
  _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
ExchangeBiasLL::computeQpResidual()
{
  if (_component == 0)
  {
    return (_eEB[_qp]*_g0[_qp]*(_alpha[_qp]*Utility::pow<2>(_mag_z[_qp])*_normals[_qp](0) + _mag_y[_qp]*(_alpha[_qp]*_mag_y[_qp]*_normals[_qp](0) - _alpha[_qp]*_mag_x[_qp]*_normals[_qp](1) + _normals[_qp](2)) - _mag_z[_qp]*(_normals[_qp](1) + _alpha[_qp]*_mag_x[_qp]*_normals[_qp](2)))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return (_eEB[_qp]*_g0[_qp]*(_alpha[_qp]*Utility::pow<2>(_mag_z[_qp])*_normals[_qp](1) + _mag_x[_qp]*(-(_alpha[_qp]*_mag_y[_qp]*_normals[_qp](0)) + _alpha[_qp]*_mag_x[_qp]*_normals[_qp](1) - _normals[_qp](2)) + _mag_z[_qp]*(_normals[_qp](0) - _alpha[_qp]*_mag_y[_qp]*_normals[_qp](2)))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 2)
  {
    return (_eEB[_qp]*_g0[_qp]*(-(_mag_y[_qp]*(_normals[_qp](0) + _alpha[_qp]*_mag_z[_qp]*_normals[_qp](1))) + _alpha[_qp]*Utility::pow<2>(_mag_y[_qp])*_normals[_qp](2) + _mag_x[_qp]*(-(_alpha[_qp]*_mag_z[_qp]*_normals[_qp](0)) + _normals[_qp](1) + _alpha[_qp]*_mag_x[_qp]*_normals[_qp](2)))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else
    return 0.0;
}

Real
ExchangeBiasLL::computeQpJacobian()
{
  if (_component == 0)
  {
    return -((_alpha[_qp]*_eEB[_qp]*_g0[_qp]*(_mag_y[_qp]*_normals[_qp](1) + _mag_z[_qp]*_normals[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]));
  }
  else if (_component == 1)
  {
    return -((_alpha[_qp]*_eEB[_qp]*_g0[_qp]*(_mag_x[_qp]*_normals[_qp](0) + _mag_z[_qp]*_normals[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]));
  }
  else if (_component == 2)
  {
    return -((_alpha[_qp]*_eEB[_qp]*_g0[_qp]*(_mag_x[_qp]*_normals[_qp](0) + _mag_y[_qp]*_normals[_qp](1))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]));
  }
  else
    return 0.0;
}

Real
ExchangeBiasLL::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (_eEB[_qp]*_g0[_qp]*(2*_alpha[_qp]*_mag_y[_qp]*_normals[_qp](0) - _alpha[_qp]*_mag_x[_qp]*_normals[_qp](1) + _normals[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
      return (_eEB[_qp]*_g0[_qp]*(2*_alpha[_qp]*_mag_z[_qp]*_normals[_qp](0)*_phi[_j][_qp] - (_normals[_qp](1) + _alpha[_qp]*_mag_x[_qp]*_normals[_qp](2))*_phi[_j][_qp])*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _mag_x_var)
    {
      return -((_eEB[_qp]*_g0[_qp]*(_alpha[_qp]*_mag_y[_qp]*_normals[_qp](0) - 2*_alpha[_qp]*_mag_x[_qp]*_normals[_qp](1) + _normals[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return (_eEB[_qp]*_g0[_qp]*(_normals[_qp](0) + 2*_alpha[_qp]*_mag_z[_qp]*_normals[_qp](1) - _alpha[_qp]*_mag_y[_qp]*_normals[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _mag_x_var)
    {
      return (_eEB[_qp]*_g0[_qp]*(-(_alpha[_qp]*_mag_z[_qp]*_normals[_qp](0)) + _normals[_qp](1) + 2*_alpha[_qp]*_mag_x[_qp]*_normals[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag_y_var)
    {
      return -((_eEB[_qp]*_g0[_qp]*(_normals[_qp](0) + _alpha[_qp]*_mag_z[_qp]*_normals[_qp](1) - 2*_alpha[_qp]*_mag_y[_qp]*_normals[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
