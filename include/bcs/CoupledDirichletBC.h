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

#ifndef COUPLEDDIRICHLETBC_H
#define COUPLEDDIRICHLETBC_H

#include "NodalBC.h"

/**
 * Implements a coupled Dirichlet BC where u = coupled_var on the boundary.
 */
class CoupledDirichletBC : public NodalBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  CoupledDirichletBC(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

private:
  /**
   * Holds the values at the quadrature points
   * of a coupled variable.
   */
  const VariableValue & _coupled_var_val;
};

#endif // COUPLEDDIRICHLETBC_H
