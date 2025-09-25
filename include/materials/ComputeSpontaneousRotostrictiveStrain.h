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

#ifndef COMPUTESPONTANEOUSROTOSTRICTIVESTRAIN_H
#define COMPUTESPONTANEOUSROTOSTRICTIVESTRAIN_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "ComputeSpontaneousRotostrictiveStrain.h"
#include "ComputeEigenstrainBase.h"

/**
 * ComputeSpontaneousRotostrictiveStrain the base class for computing spontaneous rotostrictive strain contributions (cubic)
 */
class ComputeSpontaneousRotostrictiveStrain : public ComputeEigenstrainBase
{
public:
    ComputeSpontaneousRotostrictiveStrain(const InputParameters & parameters);

  static InputParameters validParams();
  void computeQpEigenstrain();
private:
  const VariableValue & _antiphase_A_x;
  const VariableValue & _antiphase_A_y;
  const VariableValue & _antiphase_A_z;


  const MaterialProperty<Real> & _R11;
  const MaterialProperty<Real> & _R12;
  const MaterialProperty<Real> & _R44;
  std::vector<Real> _vals;
  RankTwoTensor _rotostrictive_strain;
};

#endif //COMPUTESPONTANEOUSROTOSTRICTIVESTRAIN_H
