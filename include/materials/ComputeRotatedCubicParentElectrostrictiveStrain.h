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

#ifndef COMPUTEROTATEDCUBICPARENTELECTROSTRICTIVESTRAIN_H
#define COMPUTEROTATEDCUBICPARENTELECTROSTRICTIVESTRAIN_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "ComputeEigenstrainBase.h"

/**
 * ComputeRotatedCubicParentElectrostrictiveStrain the base class for computing spontaneous polar strain contributions (cubic)
 */
class ComputeRotatedCubicParentElectrostrictiveStrain : public ComputeEigenstrainBase
{
public:
    ComputeRotatedCubicParentElectrostrictiveStrain(const InputParameters & parameters);

  static InputParameters validParams();
  void computeQpEigenstrain();
private:
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _Q11;
  const MaterialProperty<Real> & _Q12;
  const MaterialProperty<Real> & _Q44;
  const MaterialProperty<Real> & _twist;
  std::vector<Real> _vals;
  RankTwoTensor _rotated_cubic_parent_electrostrictive_strain;
};

#endif //ROTATEDCOMPUTECUBICPARENTELECTROSTRICTIVESTRAIN_H
