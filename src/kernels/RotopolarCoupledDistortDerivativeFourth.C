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

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "RotopolarCoupledDistortDerivativeFourth.h"

template<>
InputParameters validParams<RotopolarCoupledDistortDerivativeFourth>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt");
  params.addCoupledVar("antiferrodis_A_y", 0.0, "The y component of the antiferrodistortive tilt");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt");
  params.addRequiredParam<Real>("t11", "The coefficients of the coupled-Landau expansion");
  params.addRequiredParam<Real>("t12", "The coefficients of the coupled-Landau expansion");
  params.addRequiredParam<Real>("t44", "The coefficients of the coupled-Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotopolarCoupledDistortDerivativeFourth::RotopolarCoupledDistortDerivativeFourth(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _t11(getParam<Real>("t11")),
   _t12(getParam<Real>("t12")),
   _t44(getParam<Real>("t44")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotopolarCoupledDistortDerivativeFourth::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (-2*_antiferrodis_A_x[_qp]*std::pow(_polar_x[_qp],2)*_t11 - 2*_antiferrodis_A_x[_qp]*(std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))*_t12 - (_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp])*_t44);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (-2*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],2)*_t11 - 2*_antiferrodis_A_y[_qp]*(std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))*_t12 - (_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t44);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (-2*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],2)*_t11 - 2*_antiferrodis_A_z[_qp]*(std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*_t12 - (_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t44);
  }
  else
    return 0.0;
}

Real
RotopolarCoupledDistortDerivativeFourth::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*std::pow(_polar_x[_qp],2)*_t11 - 2*(std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))*_t12);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*std::pow(_polar_y[_qp],2)*_t11 - 2*(std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))*_t12);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*std::pow(_polar_z[_qp],2)*_t11 - 2*(std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*_t12);
  }
  else
    return 0.0;
}

Real
RotopolarCoupledDistortDerivativeFourth::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_x_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_x[_qp]*_polar_x[_qp]*_t11 - (_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t44);
      }
    else if (jvar == _polar_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_x[_qp]*_polar_y[_qp]*_t12 - _antiferrodis_A_y[_qp]*_polar_x[_qp]*_t44);
      }
    else if (jvar == _polar_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_x[_qp]*_polar_z[_qp]*_t12 - _antiferrodis_A_z[_qp]*_polar_x[_qp]*_t44);
      }
    else if (jvar == _antiferrodis_A_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-(_polar_x[_qp]*_polar_y[_qp]*_t44));
      }
    else if (jvar == _antiferrodis_A_z_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-(_polar_x[_qp]*_polar_z[_qp]*_t44));
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
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_t12 - _antiferrodis_A_x[_qp]*_polar_y[_qp]*_t44);
      }
    else if (jvar == _polar_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_t11 - (_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t44);
      }
    else if (jvar == _polar_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_y[_qp]*_polar_z[_qp]*_t12 - _antiferrodis_A_z[_qp]*_polar_y[_qp]*_t44);
      }
    else if (jvar == _antiferrodis_A_x_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-(_polar_x[_qp]*_polar_y[_qp]*_t44));
      }
    else if (jvar == _antiferrodis_A_z_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-(_polar_y[_qp]*_polar_z[_qp]*_t44));
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
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t12 - _antiferrodis_A_x[_qp]*_polar_z[_qp]*_t44);
      }
    else if (jvar == _polar_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t12 - _antiferrodis_A_y[_qp]*_polar_z[_qp]*_t44);
      }
    else if (jvar == _polar_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-4*_antiferrodis_A_z[_qp]*_polar_z[_qp]*_t11 - (_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp])*_t44);
      }
    else if (jvar == _antiferrodis_A_x_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-(_polar_x[_qp]*_polar_z[_qp]*_t44));
      }
    else if (jvar == _antiferrodis_A_y_var)
      {
        return _test[_i][_qp] * _phi[_j][_qp] * (-(_polar_y[_qp]*_polar_z[_qp]*_t44));
      }
    else
      {
        return 0.0;
      }
  }
  else
    return 0.0;
}
