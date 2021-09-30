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

   You should have received a co_antiferrodis_A_y[_i] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/
#pragma once

#include "MultiAppTransfer.h"

// Forward declarations
class MultiAppFieldAddTransfer;
class MooseVariableFieldBase;
namespace libMesh
{
class DofObject;
}

template <>
InputParameters validParams<MultiAppFieldAddTransfer>();

/**
 *  intermediary class that allows variable names as inputs
 */
class MultiAppFieldAddTransfer : public MultiAppTransfer
{
public:
  static InputParameters validParams();

  MultiAppFieldAddTransfer(const InputParameters & parameters);

  virtual void initialSetup();

protected:
  /**
   * Performs the transfer of a variable between two problems if they have the same mesh.
   */
  void transfer(FEProblemBase & to_problem, FEProblemBase & from_problem);

  /**
   * Performs the transfer of values between a node or element.
   */
  void transferDofObject(libMesh::DofObject * to_object,
                         libMesh::DofObject * from_object,
                         MooseVariableFieldBase & to_var,
                         MooseVariableFieldBase & from_var,
                         NumericVector<Number> & to_solution,
                         NumericVector<Number> & from_solution);

  /// Virtual function defining variables to be transferred
  virtual std::vector<VariableName> getFromVarNames() const = 0;
  /// Virtual function defining variables to transfer to
  virtual std::vector<AuxVariableName> getToVarNames() const = 0;
};
