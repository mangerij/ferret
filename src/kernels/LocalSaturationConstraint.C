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

#include "LocalSaturationConstraint.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", LocalSaturationConstraint);

template<>
InputParameters validParams<LocalSaturationConstraint>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetic vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetic vector");
  params.addRequiredParam<Real>("var_mag", "Small perturbation of magnetic constraint to eliminate degeneracy");
  params.addRequiredParam<Real>("sat_penalty", "the residual penalty for not obeying local saturation");
  return params;
}

LocalSaturationConstraint::LocalSaturationConstraint(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _var_mag(getParam<Real>("var_mag")),
  _sat_penalty(getParam<Real>("sat_penalty"))
{
}

Real
LocalSaturationConstraint::computeQpResidual()
{
  if (_component == 0)
  {
    return _sat_penalty * _test[_i][_qp] * ( Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]) - 1.0);
  }
  else if (_component == 1)
  {
    return _sat_penalty * _test[_i][_qp] * ( Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]) - 1.0);
  }
  else if (_component == 2)
  {
    return _sat_penalty * _test[_i][_qp] * ( Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]) - 1.0);
  }
  else
    return 0.0;
}

Real
LocalSaturationConstraint::computeQpJacobian()
{
  if (_component == 0)
  {
    return _sat_penalty * _test[_i][_qp] * _phi[_j][_qp] * (2.0 * _mag_x[_qp]);
  }
  else if (_component == 1)
  {
    return _sat_penalty * _test[_i][_qp] * _phi[_j][_qp] * (2.0 * _mag_y[_qp]);
  }
  else if (_component == 2)
  {
    return _sat_penalty * _test[_i][_qp] * _phi[_j][_qp] * (2.0 * _mag_z[_qp]);
  }
  else
    return 0.0;
}

//Real
//LocalSaturationConstraint::computeQpOffDiagJacobian(unsigned int jvar)
//{
//  return 0.0;
//}
