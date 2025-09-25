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

//* This file is adapted from part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MultiAppFieldAddTransfer.h"

/**
 * Add the value to the target domain from the nearest node in the source domain.
 */
class MultiAppAddTransfer : public MultiAppFieldAddTransfer
{
public:
  static InputParameters validParams();

  MultiAppAddTransfer(const InputParameters & parameters);

  /**
   * Performs the transfer of a variable (Nonlinear or Auxiliary) to/from the Multiapp.
   */
  virtual void execute() override;

protected:
  virtual std::vector<VariableName> getFromVarNames() const override { return _from_var_names; }
  virtual std::vector<AuxVariableName> getToVarNames() const override { return _to_var_names; }

  /// Name of variables transfering from
  const std::vector<VariableName> _from_var_names;
  /// Name of variables transfering to
  const std::vector<AuxVariableName> _to_var_names;

  /// This values are used if a derived class only supports one variable
  VariableName _from_var_name;
  AuxVariableName _to_var_name;
};
