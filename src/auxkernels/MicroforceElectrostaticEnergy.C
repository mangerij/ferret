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

#include "MicroforceElectrostaticEnergy.h"

registerMooseObject("FerretApp", MicroforceElectrostaticEnergy);

InputParameters MicroforceElectrostaticEnergy::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Computes the microforce due to the local electrostatic coupling.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization vector");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization vector");
  params.addCoupledVar("polar_z", "The z component of the polarization vector");
  params.addRequiredCoupledVar("potential_E_int", "The internal electrostatic potential variable");
  params.addCoupledVar("potential_E_ext", 0.0, "The external electrostatic potential variable");
  return params;
}

MicroforceElectrostaticEnergy::MicroforceElectrostaticEnergy(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _potential_E_int_grad(coupledGradient("potential_E_int")),
   _potential_E_ext_grad(coupledGradient("potential_E_ext"))
{
}

Real
MicroforceElectrostaticEnergy::computeValue()
{
//factor of 2? sign? units?
  if (_component == 0)
  {
   return (_potential_E_int_grad[_qp](_component) + _potential_E_ext_grad[_qp](_component));
  }
  else if (_component == 1)
  {
    return (_potential_E_int_grad[_qp](_component) + _potential_E_ext_grad[_qp](_component));
  }
  else if (_component == 2)
  {
    return (_potential_E_int_grad[_qp](_component) + _potential_E_ext_grad[_qp](_component));
  }
  else
    return 0.0;
}
