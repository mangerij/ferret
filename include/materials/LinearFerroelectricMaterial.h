/**
 * @file   LinearFerroelectricMaterial.h
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Sun Nov 24 22:54:31 2013
 * @brief LinearFerroelectricMaterial consider linear relation between linear elastic Tensor and electrostrictive tensor
 *       i.e. q_ijkl=2*C_ijmn * Q_mnkl
 */
#ifndef LINEARFERROELECTRICMATERIAL_H
#define LINEARFERROELECTRICMATERIAL_H

#include "TensorMechanicsMaterial.h"
#include "ElectrostrictiveTensorR4.h"
#include "libmesh/quadrature.h"
//#include "RotationTensor.h"

//Forward Declarations
class LinearFerroelectricMaterial;

template<>
InputParameters validParams<LinearFerroelectricMaterial>();

class LinearFerroelectricMaterial : public TensorMechanicsMaterial
{
public:
  LinearFerroelectricMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpStrain();
  virtual void computeQpStress();
  virtual RankTwoTensor computeStressFreeStrain();
  virtual void computeQpElectrostrictiveCoefficients();
  virtual void computeProperties();
//  virtual void computeValue();

  MaterialProperty<RankFourTensor> & _electrostrictivecoefficients;
  MaterialProperty<ElectrostrictiveTensorR4> & _electrostrictive_tensor;

  /// determines the translation from C_ijkl to the Rank-4 tensor
  RankFourTensor::FillMethod _fill_method;

  // vectors to get the input values
  std::vector<Real> _Qmnkl_vector;

//  RealVectorValue _Electrostrictive_Euler_angles;

  /// Individual material information
  RankFourTensor _Qmnkl; //electrostrictive coefficientsw
  ElectrostrictiveTensorR4 _qijkl; // q_ijkl = 2 * C_ijmn * Q_mnkl

private:
  const VariableValue & _T;

  const Real _T0;
  Real _thermal_expansion_coeff;

  std::vector<Real> _applied_strain_vector;
  RankTwoTensor _applied_strain_tensor;
};


#endif //LINEARFERROELECTRICMATERIAL_H
