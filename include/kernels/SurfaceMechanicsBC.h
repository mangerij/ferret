/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef SURFACEMECHANICSBC_H
#define SURFACEMECHANICSBC_H

#include "IntegratedBC.h"

//Forward Declarations
class SurfaceMechanicsBC;
class RankTwoTensor;
class RankFourTensor;

template<>
InputParameters validParams<SurfaceMechanicsBC>();

class SurfaceMechanicsBC : public IntegratedBC
{
public:

  SurfaceMechanicsBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual RankTwoTensor computeQpProjection();

  virtual RankTwoTensor computeQpRankTwoRotation(unsigned int tensor_component);

  virtual RankFourTensor computeQpRankFourRotation(unsigned int tensor_component);

  const unsigned int _component;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;

  std::vector<Real> _Csijkl_vector; //surfaceFillFromInputVector method only takes std::vector<Real>. Might need to change..

  RealVectorValue _S_k_vector;

  Real _taus;
  
  RealVectorValue _surface_euler_angles;

};
#endif // SURFACEMECHANICSBC_H
