/* AUTO-DERIVED from include/postprocessors/CubicParentElasticEnergy.h by substituting the polar order parameter for the
 * antiphase (AFD) tilt and the electrostrictive tensor Q for the rotostrictive
 * tensor R.  Same free energy, 1/2 (eps - eps0):C:(eps - eps0), with
 * eps0 = R.AA instead of Q.PP.  Generated 2026-09-14. */
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

#ifndef CUBICPARENTELASTICAENERGY_H
#define CUBICPARENTELASTICAENERGY_H

#include "ElementIntegralPostprocessor.h"

class CubicParentElasticAEnergy : public ElementIntegralPostprocessor
{
public:
  CubicParentElasticAEnergy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpIntegral();

private:
  const VariableValue & _antiphase_A_x;
  const VariableValue & _antiphase_A_y;
  const VariableValue & _antiphase_A_z;
  const MaterialProperty<Real> & _R11;
  const MaterialProperty<Real> & _R12;
  const MaterialProperty<Real> & _R44;
  const MaterialProperty<Real> & _C11;
  const MaterialProperty<Real> & _C12;
  const MaterialProperty<Real> & _C44;
  const Real _energy_scale;
  const std::string _base_name;
  const MaterialProperty<RankTwoTensor> & _strain;
};

#endif
