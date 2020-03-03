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

#include "HigherOrderMagConstraintCartLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", HigherOrderMagConstraintCartLLG);

template<>
InputParameters validParams<HigherOrderMagConstraintCartLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the constraint |M|-Ms = 0 using the method of lagrange multipliers");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the magnetic vector");
  params.addRequiredCoupledVar("lambda", "The lagrange multiplier");
  return params;
}

HigherOrderMagConstraintCartLLG::HigherOrderMagConstraintCartLLG(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _lambda_var(coupled("lambda")),
   _lambda(coupledValue("lambda"))
{
}


Real
HigherOrderMagConstraintCartLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return 0.0;
  }
    else if (_component == 1)
  {
    return 0.0;
  }
    else if (_component == 2)
  {
    return 0.0;
  }
  else
    return 0.0;
}

Real
HigherOrderMagConstraintCartLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return 0.0;
  }
    else if (_component == 1)
  {
    return 0.0;
  }
    else if (_component == 2)
  {
    return 0.0;
  }
    else
      return 0.0;
}

Real
HigherOrderMagConstraintCartLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return 0.0;
    }
    else if (jvar == _mag_z_var)
    {
      return 0.0;
    }
    else if (jvar == _lambda_var)
    {
      return 0.0;
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
      return 0.0;
    }
    else if (jvar == _mag_z_var)
    {
      return 0.0;
    }
    else if (jvar == _lambda_var)
    {
      return 0.0;
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
      return 0.0;
    }
    else if (jvar == _mag_y_var)
    {
      return 0.0;
    }
    else if (jvar == _lambda_var)
    {
      return 0.0;
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
