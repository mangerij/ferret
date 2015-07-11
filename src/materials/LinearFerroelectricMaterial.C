/***
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
  // params.addParam<Real>("electrostrictive_euler_angle_1", 0.0, "Euler angle in direction 1 for electrostrictive tensor");
  // params.addParam<Real>("electrostrictive_euler_angle_2", 0.0, "Euler angle in direction 2 for electrostrictive tensor");
  // params.addParam<Real>("electrostrictive_euler_angle_3", 0.0, "Euler angle in direction 3 for electrostrictive tensor");
  return params;
}

LinearFerroelectricMaterial::LinearFerroelectricMaterial(const std::string & name,
                                 InputParameters parameters) :

  LinearElasticMaterial(name, parameters),
  _electrostrictivecoefficients(declareProperty<ElasticityTensorR4>("electrostrictivecoefficients")),
  _electrostrictive_tensor(declareProperty<ElectrostrictiveTensorR4>("electrostrictive_tensor")),
  _fill_method((RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
  _Qmnkl_vector(getParam<std::vector<Real> >("Q_mnkl")),
//  _Electrostrictive_Euler_angles(getParam<Real>("electrostrictive_euler_angle_1"),
//                                 getParam<Real>("electrostrictive_euler_angle_2"),
//                                 getParam<Real>("electrostrictive_euler_angle_3")),
  _Qmnkl(),
  _qijkl()
{
  _Qmnkl.fillFromInputVector(_Qmnkl_vector,_fill_method);
//_qijkl.computeValue(_Cijkl,_Qmnkl);
}

void
LinearFerroelectricMaterial::computeQpElectrostrictiveCoefficients()
{
  // eR type: RealTensorValue
//  RotationTensor eR(_Electrostrictive_Euler_angles);

  _electrostrictivecoefficients[_qp] = _Qmnkl;
//  _electrostrictivecoefficients[_qp].rotate(eR) //construct rotated Q_mnkl

}

void
LinearFerroelectricMaterial::computeProperties()
{
  {
  // Moose::out << "\n Performing C_ijmn Q_mnkl contraction on all the quadrature points?";
  for(_qp = 0; _qp < _qrule->n_points(); ++_qp)
  _electrostrictive_tensor[_qp].computeProduct(_Cijkl,_Qmnkl);
  }
}

// note that the elasticity properties are loaded into the input file with
// a different materials block
