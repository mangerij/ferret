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

#include "AnisotropyEnergy.h"
#include "libmesh/utility.h"

class AnisotropyEnergy;

template<>
InputParameters validParams<AnisotropyEnergy>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Adds a residual contribution for an arbitrary anisotropy quadratic in the polarization field.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("K", 1.0, "the anisotropy energy");
  return params;
}

AnisotropyEnergy::AnisotropyEnergy(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _K(getParam<Real>("K"))
{
}

Real
AnisotropyEnergy::computeQpResidual()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  return 2.0 * _K * w(_component) * _test[_i][_qp] * Utility::pow<3>(_len_scale);
}

Real
AnisotropyEnergy::computeQpJacobian()
{
  return 2.0 * _K * _phi[_j][_qp] * _test[_i][_qp] * Utility::pow<3>(_len_scale);
}
