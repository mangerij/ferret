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

#include "DzyaloshinskiiMagDerivative.h"

registerMooseObject("FerretApp", DzyaloshinskiiMagDerivative);

template<>
InputParameters validParams<DzyaloshinskiiMagDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution in magnetoelectric coupling due to antiferrodistortive and magnetic ordering.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the antiferromagnetic vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive vector");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive vector");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive vector");
  params.addRequiredParam<Real>("chiP", "The magnetic susceptibility");
  params.addRequiredParam<Real>("hD", "The Dzyaloshinkii field");
  params.addRequiredParam<Real>("M0", "The scaling of the magnetic vector");
  params.addRequiredParam<Real>("A0", "The scaling of the antiferrodistortive tilt vector");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

DzyaloshinskiiMagDerivative::DzyaloshinskiiMagDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
  _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
  _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
  _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
  _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
  _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
  _chiP(getParam<Real>("chiP")),
  _hD(getParam<Real>("hD")),
  _M0(getParam<Real>("M0")),
  _A0(getParam<Real>("A0")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
DzyaloshinskiiMagDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _chiP * std::pow(_hD,2.0)*(_antiferrodis_A_z[_qp]*_mag_y[_qp] - _antiferrodis_A_y[_qp]*_mag_z[_qp])*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]);
  }
  else if (_component == 1)
  {
    return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _chiP * std::pow(_hD,2.0)*(-(_antiferrodis_A_z[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_z[_qp])*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]);
  }
  else if (_component == 2)
  {
    return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _chiP * std::pow(_hD,2.0)*(_antiferrodis_A_y[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_y[_qp])*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]);
  }
  else
    return 0.0;
}

Real
DzyaloshinskiiMagDerivative::computeQpJacobian()
{
  //now since this is a double cross product, on-diagonal terms appear.
  if (_component == 0)
  {
    return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_x[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_z[_qp]*_mag_y[_qp] - _antiferrodis_A_y[_qp]*_mag_z[_qp]));
  }
  else if (_component == 1)
  {
    return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_y[_qp]*_chiP*std::pow(_hD,2)*(-(_antiferrodis_A_z[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_z[_qp]));
  }
  else if (_component == 2)
  {
    return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_z[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_y[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_y[_qp]));
  }
  else
    return 0.0;
}

Real
DzyaloshinskiiMagDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_y[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_z[_qp]*_mag_y[_qp] - _antiferrodis_A_y[_qp]*_mag_z[_qp]) + _antiferrodis_A_z[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_z[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_z[_qp]*_mag_y[_qp] - _antiferrodis_A_y[_qp]*_mag_z[_qp]) - _antiferrodis_A_y[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_x[_qp]*(_antiferrodis_A_z[_qp]*_mag_y[_qp] - _antiferrodis_A_y[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_y[_qp]*(_antiferrodis_A_z[_qp]*_mag_y[_qp] - _antiferrodis_A_y[_qp]*_mag_z[_qp]) - _chiP*std::pow(_hD,2)*_mag_z[_qp]*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_z[_qp]*(_antiferrodis_A_z[_qp]*_mag_y[_qp] - _antiferrodis_A_y[_qp]*_mag_z[_qp]) + _chiP*std::pow(_hD,2)*_mag_y[_qp]*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
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
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_x[_qp]*_chiP*std::pow(_hD,2)*(-(_antiferrodis_A_z[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_z[_qp]) - _antiferrodis_A_z[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_z[_qp]*_chiP*std::pow(_hD,2)*(-(_antiferrodis_A_z[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_z[_qp]) + _antiferrodis_A_x[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_x[_qp]*(-(_antiferrodis_A_z[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_z[_qp]) + _chiP*std::pow(_hD,2)*_mag_z[_qp]*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_y[_qp]*(-(_antiferrodis_A_z[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_z[_qp]*(-(_antiferrodis_A_z[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_z[_qp]) - _chiP*std::pow(_hD,2)*_mag_x[_qp]*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
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
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_x[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_y[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_y[_qp]) + _antiferrodis_A_y[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_antiferrodis_A_y[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_y[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_y[_qp]) - _antiferrodis_A_x[_qp]*_chiP*std::pow(_hD,2)*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_x[_qp]*(_antiferrodis_A_y[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_y[_qp]) - _chiP*std::pow(_hD,2)*_mag_y[_qp]*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*_mag_y[_qp]*(_antiferrodis_A_y[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_y[_qp]) + _chiP*std::pow(_hD,2)*_mag_x[_qp]*(_antiferrodis_A_x[_qp]*_mag_x[_qp] + _antiferrodis_A_y[_qp]*_mag_y[_qp] + _antiferrodis_A_z[_qp]*_mag_z[_qp]));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return (1.0 / (4 * _M0 * _M0 * _A0 * _A0) ) * _test[_i][_qp] * _phi[_j][_qp] * (_chiP*std::pow(_hD,2)*(_antiferrodis_A_y[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_y[_qp])*_mag_z[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
