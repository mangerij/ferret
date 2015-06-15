/**
 * @file   LinearFerroelectricMaterial.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15 2015
 * @brief  LinearFerroelectricMaterial consider linear relation between linear elastic Tensor and electrostrictive tensor
 *       i.e. q_ijkl=2*C_ijmm * Q_mnkl
 */


#include "LinearFerroelectricMaterial.h"

template<>
InputParameters validParams<LinearFerroelectricMaterial>()
{
  InputParameters params = validParams<LinearElasticMaterial>();
  params.addRequiredParam<std::vector<Real> >("Q_mnkl", "electrostrictive coefficients(vector)");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

LinearFerroelectricMaterial::LinearFerroelectricMaterial(const std::string & name,
                                 InputParameters parameters) :

    LinearElasticMaterial(name, parameters),
    _electrostrictivecoefficients(declareProperty<ElasticityTensorR4>("electrostrictivecoefficients")),
    _electrostrictive_tensor(declareProperty<ElectrostrictiveTensorR4>("electrostrictive_tensor")),
    _fill_method((RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _Qmnkl_vector(getParam<std::vector<Real> >("Q_mnkl")),
    _Qmnkl(),
    _qijkl()
{
  _Qmnkl.fillFromInputVector(_Qmnkl_vector,_fill_method);
//_qijkl.computeValue(_Cijkl,_Qmnkl);
}

void
LinearFerroelectricMaterial::computeQpElectrostrictiveCoefficients()
{
  _electrostrictivecoefficients[_qp] = _Qmnkl;
}

void
LinearFerroelectricMaterial::computeProperties()
{
  for(_qp=0; _qp < _qrule->n_points(); ++_qp)
  _electrostrictive_tensor[_qp].computeProduct(_Cijkl,_Qmnkl);
}
//  for(_qp=0; _qp < _qrule->n_points(); ++_qp) used to find a solution?

// void
// LinearFerroelectricMaterial::computeProperties()
// {
//  for(_qp=0; _qp < _qrule->n_points(); ++_qp)
//  {
//    _electrostrictive_tensor[_qp]=_qijkl;
//  }
//  LinearElasticMaterial::computeProperties();
// }
