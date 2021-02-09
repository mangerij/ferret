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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "CoupledDirichletBC.h"

registerMooseObject("FerretApp", CoupledDirichletBC);

template <>
InputParameters
validParams<CoupledDirichletBC>()
{
  InputParameters params = validParams<NodalBC>();

  params.addRequiredCoupledVar("coupled_var", "Value on the Boundary");
  return params;
}

CoupledDirichletBC::CoupledDirichletBC(const InputParameters & parameters)
  : NodalBC(parameters),

    /**
     * Get a reference to the coupled variable's values.
     */
    _coupled_var_val(coupledValue("coupled_var"))
{
}

Real
CoupledDirichletBC::computeQpResidual()
{
  return _u[_qp] - _coupled_var_val[_qp];
}
