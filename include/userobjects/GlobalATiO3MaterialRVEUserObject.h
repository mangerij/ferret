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

#ifndef GLOBALATIO3MATERIALRVEUSEROBJECT_H
#define GLOBALATIO3MATERIALRVEUSEROBJECT_H

#include "ElementUserObject.h"
#include "GlobalStrainUserObjectInterface.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class GlobalATiO3MaterialRVEUserObject : public ElementUserObject, public GlobalStrainUserObjectInterface
{
public:
  GlobalATiO3MaterialRVEUserObject(const InputParameters & parameters);

  static InputParameters validParams();

  void initialize() override;
  void execute() override;
  void threadJoin(const UserObject & uo) override;
  void finalize() override;
  virtual const RankTwoTensor & getResidual() const override;
  virtual const RankFourTensor & getJacobian() const override;
  virtual const VectorValue<bool> & getPeriodicDirections() const override;

  /**
   * Calculate additional applied stresses
   */
  virtual void computeAdditionalStress(){};

protected:
  std::string _base_name;

  const MaterialProperty<RankFourTensor> & _dstress_dstrain;
  const MaterialProperty<RankTwoTensor> & _stress;

  RankTwoTensor _applied_stress_tensor;
  RankTwoTensor _residual;
  RankFourTensor _jacobian;

  const unsigned int _dim;
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  VectorValue<bool> _periodic_dir;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _C11;
  const MaterialProperty<Real> & _C12;
  const MaterialProperty<Real> & _C44;
  const MaterialProperty<Real> & _Q11;
  const MaterialProperty<Real> & _Q12;
  const MaterialProperty<Real> & _Q44;

};

#endif // GLOBALATIO3MATERIALRVEUSEROBJECT_H
