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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "FourierHeat.h"
#include "Material.h"

registerMooseObject("FerretApp", FourierHeat);

InputParameters
FourierHeat::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to k*deltaT = 0");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<MaterialPropertyName>("thC", "Thermal Conductivity", "Property name of the thermal conductivity material property");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

FourierHeat::FourierHeat(const InputParameters & parameters)
  :AuxKernel(parameters),
   _thC(getMaterialProperty<Real>("thC")),
   _component(getParam<unsigned int>("component")),
   _T_var(coupled("T")),
   _T(coupledValue("T")),
   _T_grad(coupledGradient("T")),
   _len_scale(getParam<Real>("len_scale"))
{
}
Real
FourierHeat::computeValue()
{
  return -_thC[_qp] * _T_grad[_qp](_component) * _len_scale;
}
