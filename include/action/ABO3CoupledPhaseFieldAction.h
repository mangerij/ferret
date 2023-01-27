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

#pragma once

#include "Action.h"
#include "libmesh/fe_type.h"

class ABO3CoupledPhaseFieldAction : public Action
{
public:
  ABO3CoupledPhaseFieldAction(const InputParameters & params);
  static InputParameters validParams();
  virtual void act() override;

protected:
  const bool _coupled_problem;
  const bool _polar_time_dependence;
  const bool _u_time_dependence;
  const bool _phi_time_dependence;
  const bool _is_renormalized;
  const bool _is_permittivity_anisotropic;
};
