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

#pragma once

#include "Action.h"
#include "libmesh/fe_type.h"

class CubicParentFEPhaseFieldAction : public Action
{
public:
  CubicParentFEPhaseFieldAction(const InputParameters & params);
  static InputParameters validParams();
  virtual void act() override;

protected:
  void addConstantMaterial(const std::string & prop_param,
                           const std::string & val_param,
                           const std::string & object_name);

  void setPolarCoupling(InputParameters & params) const;

  const std::vector<NonlinearVariableName> _polar_vars;
  const NonlinearVariableName _potential_var;
  const bool _electrostatics;
  const bool _elastic;
  const bool _polar_time_dependence;
  const bool _add_variables;
  const bool _add_elastic_materials;
  const bool _add_elastic_energy_aux;
  const bool _add_postprocessors;
  const std::string _eigenstrain_name;
};
