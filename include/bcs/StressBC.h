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

#ifndef STRESSBC_H
#define STRESSBC_H

#include "IntegratedBC.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class StressBC;

template<>
InputParameters validParams<StressBC>();

class StressBC : public IntegratedBC
{
public:
  StressBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  std::vector<Real> _stress_vector;

  RankTwoTensor _boundary_stress;
  std::vector<const VariableValue *> _boundary_stress_vars;

  const MaterialProperty<RankFourTensor> & _Jacobian_mult;
  const int _component;
  bool _convert_to_gpa;
  Real _prefactStress;
  Real _multiplier;
};

#endif
