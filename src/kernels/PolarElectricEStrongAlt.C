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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "PolarElectricEStrongAlt.h"

class PolarElectricEStrongAlt;

template<>
InputParameters validParams<PolarElectricEStrongAlt>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

PolarElectricEStrongAlt::PolarElectricEStrongAlt(const InputParameters & parameters)
  :Kernel(parameters),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
PolarElectricEStrongAlt::computeQpResidual()
{
  Real RpolarE = 0.0;
  RpolarE += -(_polar_x_grad[_qp](0) * _test[_i][_qp] + _polar_y_grad[_qp](1) * _test[_i][_qp] + _polar_z_grad[_qp](2) * _test[_i][_qp]) * std::pow(_len_scale, 2.0); //this should be a MINUS
  ///  Moose::out << "\n R_polarE-"; std::cout << " = " << RpolarE;
  return RpolarE;
}
Real
PolarElectricEStrongAlt::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricEStrongAlt::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _polar_x_var)
    return -_grad_phi[_j][_qp](0) * _test[_i][_qp] * std::pow(_len_scale, 2.0);
  else if (jvar == _polar_y_var)
    return -_grad_phi[_j][_qp](1) * _test[_i][_qp] * std::pow(_len_scale, 2.0);
  else if (jvar == _polar_z_var)
    return -_grad_phi[_j][_qp](2) * _test[_i][_qp] * std::pow(_len_scale, 2.0);
  else
    return 0.0;
}
