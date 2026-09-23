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

#ifndef CUBICPARENTELASTICPDERIVATIVE_H
#define CUBICPARENTELASTICPDERIVATIVE_H

#include "Kernel.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

class CubicParentElasticPDerivative: public Kernel
{
public:
  CubicParentElasticPDerivative(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  /// P <-> mean strain (SCALAR global_strain) blocks, both directions. See the .C.
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

  /// d eps0 / d P_i as a tensor (unrotated cubic parent, TENSOR shear convention
  /// eps0_ij = Q44 P_i P_j, i.e. exactly what this kernel's residual differentiates).
  RankTwoTensor dEigenstrain_dP(unsigned int i) const;
  /// unit symmetric tensor for scalar mean-strain component b (MOOSE ordering
  /// xx yy zz yz xz xy, as RankTwoTensor::fillFromScalarVariable)
  static RankTwoTensor unitStrain(unsigned int b);

private:
  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _C11;
  const MaterialProperty<Real> & _C12;
  const MaterialProperty<Real> & _C44;
  const MaterialProperty<Real> & _Q11;
  const MaterialProperty<Real> & _Q12;
  const MaterialProperty<Real> & _Q44;
  const std::string _base_name;
  const MaterialProperty<RankTwoTensor> & _strain;

  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  const unsigned int _global_strain_var;
  const MaterialProperty<RankFourTensor> * _elasticity_tensor;
  DenseMatrix<Number> _ke_copy;
};
#endif //CUBICPARENTELASTICPDERIVATIVE_H
