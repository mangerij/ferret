/****************************************************************/
/* SurfaceMechanicsBC                                          */
/*     Includes surface stress contributions to the free energy.*/
/*                                                              */
/*    IntegrateBC: the independent surface elasticity constants */
/*    and Euler angles for each surface. A residual is formed   */
/*    as an integrated boundary condition, as IntegratedBC      */
/*    has access to surface integrals and surface normals       */
/*                                                              */
/****************************************************************/

#ifndef SURFACEMECHANICSBC_H
#define SURFACEMECHANICSBC_H
#include "FerretBase.h"
#include "Kernel.h"
#include "FEProblem.h"
#include "IntegratedBC.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
//#include "TensorMechanicsMaterial.h"
//LibMesh includes
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

//Forward Declarations
class SurfaceMechanicsBC;

template<>
InputParameters validParams<SurfaceMechanicsBC>();

class SurfaceMechanicsBC : public IntegratedBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  SurfaceMechanicsBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
//  virtual Real computeQpSurfaceElasticityTensor();
  virtual void computeQpProjection();
  virtual void computeQpRotation();
//  virtual void computeQpSurfaceStress();

  private:
  const unsigned int _dim;
  const unsigned int _component;


//you will want to make sure these are named according to how you want them to be
  Real _surface_euler_angle_1;
  Real _surface_euler_angle_2;
  Real _surface_euler_angle_3;
  std::vector<Real> _Csijkl_vector;
  //std::vector<Real> _tausij_vector;
  Real _taus;
  RankFourTensor _Csijkl;
  // RankTwoTensor _surface_tau;
  RealVectorValue _surface_euler_angles;
  RealVectorValue _tangent_1;
  RealVectorValue _tangent_2;

  VariableGradient & _grad_disp_x;
  VariableGradient & _grad_disp_y;
  VariableGradient & _grad_disp_z;

  RankTwoTensor _projection, _surface_strain, _surface_stress;
  RankTwoTensor _tp11, _tp22;
  RankFourTensor _t11, _t22, _t12, _t33;
  Real C0000, C1111, C0011, C0101;
  //  RankTwoTensor grad_tensor;
};

#endif
