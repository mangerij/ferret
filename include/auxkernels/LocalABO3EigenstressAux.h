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

#pragma once

#include "NodalPatchRecovery.h"

#include "RankFourTensor.h"
#include "RankTwoTensor.h"

class LocalABO3EigenstressAux;

template <> InputParameters validParams<LocalABO3EigenstressAux>();

class LocalABO3EigenstressAux : public NodalPatchRecovery {
public:
  LocalABO3EigenstressAux(const InputParameters &parameters);

protected:
  virtual Real computeValue();

  const unsigned int _i;
  const unsigned int _j;

  const VariableValue &_polar_x;
  const VariableValue &_polar_y;
  const VariableValue &_polar_z;
  const MaterialProperty<Real> &_C11;
  const MaterialProperty<Real> &_C12;
  const MaterialProperty<Real> &_C44;
  const MaterialProperty<Real> &_Q11;
  const MaterialProperty<Real> &_Q12;
  const MaterialProperty<Real> &_Q44;
};
