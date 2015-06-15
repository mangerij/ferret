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

#include "LinearElasticMaterial.h"
#include "ElectrostrictiveTensorR4.h"

//Forward Declarations
class LinearFerroelectricMaterial;

template<>
InputParameters validParams<LinearFerroelectricMaterial>();

class LinearFerroelectricMaterial : public LinearElasticMaterial
{
public:
  LinearFerroelectricMaterial(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpElectrostrictiveCoefficients();
  virtual void computeProperties();
//  virtual void computeValue();

  MaterialProperty<ElasticityTensorR4> & _electrostrictivecoefficients;
  MaterialProperty<ElectrostrictiveTensorR4> & _electrostrictive_tensor;

  /// determines the translation from C_ijkl to the Rank-4 tensor
  RankFourTensor::FillMethod _fill_method;

  // vectors to get the input values
  std::vector<Real> _Qmnkl_vector;

  /// Individual material information
  ElasticityTensorR4 _Qmnkl; //electrostrictive coefficients will set them as elasticity components right now
  ElectrostrictiveTensorR4 _qijkl; // q_ijkl=2*C_ijmn * Q_mnkl

};


#endif //LINEARFERROELECTRICMATERIAL_H
