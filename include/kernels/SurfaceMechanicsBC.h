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

/*
*    Includes surface stress contributions to the free energy
*    A residual is formed as an integrated boundary condition, as IntegratedBC
*    has access to surface integrals and surface normals
*/
class SurfaceMechanicsBC : public IntegratedBC
{
public:

  SurfaceMechanicsBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual void computeQpProjection();
  virtual void computeQpRotation();

  const unsigned int _dim;
  const unsigned int _component;

  std::vector<Real> _Csijkl_vector;
  Real _taus;
  Real _surface_Euler_angle_1;
  Real _surface_Euler_angle_2;
  Real _surface_Euler_angle_3;

  RealVectorValue _surface_Euler_angles;

  const VariableGradient & _grad_disp_x;
  const VariableGradient & _grad_disp_y;
  const VariableGradient & _grad_disp_z;
};

#endif
