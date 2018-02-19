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

#include "DampedAnisotropyLLG.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<DampedAnisotropyLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the antiferromagnetic vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the antiferromagnetic vector");
  params.addRequiredParam<Real>("alphaLL", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("Ku", "Ku");
  params.addRequiredParam<Real>("nx", "x component of the anisotropy director");
  params.addRequiredParam<Real>("ny", "y component of the anisotropy director");
  params.addRequiredParam<Real>("nz", "z component of the anisotropy director");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

DampedAnisotropyLLG::DampedAnisotropyLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _alphaLL(getParam<Real>("alphaLL")),
  _Ku(getParam<Real>("Ku")),
  _nx(getParam<Real>("nx")),
  _ny(getParam<Real>("ny")),
  _nz(getParam<Real>("nz")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
DampedAnisotropyLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_alphaLL*_Ku*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz)*(-(Utility::pow<2>(_mag_y[_qp])*_nx) + _mag_x[_qp]*_mag_y[_qp]*_ny + _mag_z[_qp]*(-(_mag_z[_qp]*_nx) + _mag_x[_qp]*_nz)));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (-2*_alphaLL*_Ku*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz)*(-(_mag_x[_qp]*_mag_y[_qp]*_nx) + Utility::pow<2>(_mag_x[_qp])*_ny + _mag_z[_qp]*(_mag_z[_qp]*_ny - _mag_y[_qp]*_nz)));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_alphaLL*_Ku*(_mag_z[_qp]*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny) - (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz));
  }
  else
    return 0.0;
}

Real
DampedAnisotropyLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_alphaLL*_Ku*(_mag_y[_qp]*_ny + _mag_z[_qp]*_nz)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz) + 2*_alphaLL*_Ku*_nx*(-(Utility::pow<2>(_mag_y[_qp])*_nx) + _mag_x[_qp]*_mag_y[_qp]*_ny + _mag_z[_qp]*(-(_mag_z[_qp]*_nx) + _mag_x[_qp]*_nz)));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*_alphaLL*_Ku*(-(_mag_x[_qp]*_nx) - _mag_z[_qp]*_nz)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz) - 2*_alphaLL*_Ku*_ny*(-(_mag_x[_qp]*_mag_y[_qp]*_nx) + Utility::pow<2>(_mag_x[_qp])*_ny + _mag_z[_qp]*(_mag_z[_qp]*_ny - _mag_y[_qp]*_nz)));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] *  (2*_alphaLL*_Ku*_nz*(_mag_z[_qp]*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny) - (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz) + 2*_alphaLL*_Ku*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz));
  }
  else
    return 0.0;
}

Real
DampedAnisotropyLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (2*_alphaLL*_Ku*(-2*_mag_y[_qp]*_nx + _mag_x[_qp]*_ny)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz) + 2*_alphaLL*_Ku*_ny*(-(Utility::pow<2>(_mag_y[_qp])*_nx) + _mag_x[_qp]*_mag_y[_qp]*_ny + _mag_z[_qp]*(-(_mag_z[_qp]*_nx) + _mag_x[_qp]*_nz)));
    }
    else if (jvar == _mag_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (2*_alphaLL*_Ku*(-2*_mag_z[_qp]*_nx + _mag_x[_qp]*_nz)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz) + 2*_alphaLL*_Ku*_nz*(-(Utility::pow<2>(_mag_y[_qp])*_nx) + _mag_x[_qp]*_mag_y[_qp]*_ny + _mag_z[_qp]*(-(_mag_z[_qp]*_nx) + _mag_x[_qp]*_nz)));
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
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_alphaLL*_Ku*(-(_mag_y[_qp]*_nx) + 2*_mag_x[_qp]*_ny)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz) - 2*_alphaLL*_Ku*_nx*(-(_mag_x[_qp]*_mag_y[_qp]*_nx) + Utility::pow<2>(_mag_x[_qp])*_ny + _mag_z[_qp]*(_mag_z[_qp]*_ny - _mag_y[_qp]*_nz)));
    }
    else if (jvar == _mag_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_alphaLL*_Ku*(2*_mag_z[_qp]*_ny - _mag_y[_qp]*_nz)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz) - 2*_alphaLL*_Ku*_nz*(-(_mag_x[_qp]*_mag_y[_qp]*_nx) + Utility::pow<2>(_mag_x[_qp])*_ny + _mag_z[_qp]*(_mag_z[_qp]*_ny - _mag_y[_qp]*_nz)));
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
      return _test[_i][_qp] * _phi[_j][_qp] * (2*_alphaLL*_Ku*_nx*(_mag_z[_qp]*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny) - (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz) + 2*_alphaLL*_Ku*(_mag_z[_qp]*_nx - 2*_mag_x[_qp]*_nz)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz));
    }
    else if (jvar == _mag_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (2*_alphaLL*_Ku*_ny*(_mag_z[_qp]*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny) - (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz) + 2*_alphaLL*_Ku*(_mag_z[_qp]*_ny - 2*_mag_y[_qp]*_nz)*(_mag_x[_qp]*_nx + _mag_y[_qp]*_ny + _mag_z[_qp]*_nz));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
