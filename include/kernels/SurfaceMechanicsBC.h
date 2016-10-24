/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef SURFACEMECHANICSBC_H
#define SURFACEMECHANICSBC_H

#include "FEProblem.h"
#include "IntegratedBC.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"


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
  virtual void computeQpProjection();
  virtual void computeQpRotation();

  const unsigned int _component;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;

  Real _surface_euler_angle_1;
  Real _surface_euler_angle_2;
  Real _surface_euler_angle_3;
  std::vector<Real> _Csijkl_vector;

  Real _taus;

  RankFourTensor _Csijkl;
  RealVectorValue _surface_euler_angles;

  RankTwoTensor _projection, _surface_strain, _surface_stress;
  RankTwoTensor _tp11, _tp22;
};

#endif // SURFACEMECHANICSBC_H
