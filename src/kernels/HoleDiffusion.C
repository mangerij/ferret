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


#include "HoleDiffusion.h"

registerMooseObject("FerretApp", HoleDiffusion);

template<>
InputParameters validParams<HoleDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

HoleDiffusion::HoleDiffusion(const InputParameters & parameters)
  :Kernel(parameters),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
HoleDiffusion::computeQpResidual()
{
  Real Relec = 0.0;
  Relec +=  _grad_u[_qp] * _grad_test[_i][_qp] * _len_scale;
  ///  Moose::out << "\n R_elec-"; std::cout << " = " << Relec;
  return Relec;
}

Real
HoleDiffusion::computeQpJacobian()
{
   return _grad_phi[_j][_qp] * _grad_test[_i][_qp] * _len_scale;
}
