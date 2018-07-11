//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef GLOBALATOMATERIALRVEUSEROBJECT_H
#define GLOBALATOMATERIALRVEUSEROBJECT_H

#include "ElementUserObject.h"
#include "GlobalStrainUserObjectInterface.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class GlobalATOMaterialRVEUserObject;

template <>
InputParameters validParams<GlobalATOMaterialRVEUserObject>();

class GlobalATOMaterialRVEUserObject : public ElementUserObject, public GlobalStrainUserObjectInterface
{
public:
  GlobalATOMaterialRVEUserObject(const InputParameters & parameters);

  void initialize() override;
  void execute() override;
  void threadJoin(const UserObject & uo) override;
  void finalize() override;
  virtual const RankTwoTensor & getResidual() const;
  virtual const RankFourTensor & getJacobian() const;
  virtual const VectorValue<bool> & getPeriodicDirections() const;

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
  const Real _C11;
  const Real _C12;
  const Real _C44;
  const Real _Q11;
  const Real _Q12;
  const Real _Q44;

};

#endif // GLOBALATOMATERIALRVEUSEROBJECT_H
