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

#ifndef COMPUTECUBICPARENTFLEXOELECTRICSTRAIN_H
#define COMPUTECUBICPARENTFLEXOELECTRICSTRAIN_H

#include "ComputeEigenstrainBase.h"
#include "RankTwoTensor.h"

/**
 * ComputeCubicParentFlexoelectricStrain computes the flexoelectric eigenstrain of a cubic parent phase.
 */
class ComputeCubicParentFlexoelectricStrain : public ComputeEigenstrainBase
{
public:
  ComputeCubicParentFlexoelectricStrain(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpEigenstrain() override;

private:
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
  const MaterialProperty<Real> & _s11;
  const MaterialProperty<Real> & _s12;
  const MaterialProperty<Real> & _s44;
  const MaterialProperty<Real> & _mu11;
  const MaterialProperty<Real> & _mu12;
  const MaterialProperty<Real> & _mu44;
  std::vector<Real> _vals;
  RankTwoTensor _flexoelectric_strain;
};

#endif //COMPUTECUBICPARENTFLEXOELECTRICSTRAIN_H
