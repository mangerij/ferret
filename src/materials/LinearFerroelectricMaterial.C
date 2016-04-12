/***
 * @file   LinearFerroelectricMaterial.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15 2015
 * @brief  LinearFerroelectricMaterial consider linear relation between linear elastic Tensor and electrostrictive tensor
 *       i.e. q_ijkl = 2 * C_ijmm * Q_mnkl
 */


#include "LinearFerroelectricMaterial.h"

template<>
InputParameters validParams<LinearFerroelectricMaterial>()
{
  InputParameters params = validParams<TensorMechanicsMaterial>();
  params.addRequiredParam<std::vector<Real> >("Q_mnkl", "electrostrictive coefficients(vector)");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addParam<Real>("thermal_expansion_coeff", 0, "Thermal expansion coefficient in 1/K");
  params.addParam<Real>("T0", 300, "Reference temperature for thermal expansion in K");
  params.addCoupledVar("T", 300, "Temperature in Kelvin");
  params.addParam<std::vector<Real> >("applied_strain_vector", "Applied strain: e11, e22, e33, e23, e13, e12");
  // params.addParam<Real>("electrostrictive_euler_angle_1", 0.0, "Euler angle in direction 1 for electrostrictive tensor");
  // params.addParam<Real>("electrostrictive_euler_angle_2", 0.0, "Euler angle in direction 2 for electrostrictive tensor");
  // params.addParam<Real>("electrostrictive_euler_angle_3", 0.0, "Euler angle in direction 3 for electrostrictive tensor");
  return params;
}

LinearFerroelectricMaterial::LinearFerroelectricMaterial(const InputParameters & parameters) :

  TensorMechanicsMaterial(parameters),
  _electrostrictivecoefficients(declareProperty<RankFourTensor>("electrostrictivecoefficients")),
  _electrostrictive_tensor(declareProperty<ElectrostrictiveTensorR4>("electrostrictive_tensor")),
  _fill_method((RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
  _Qmnkl_vector(getParam<std::vector<Real> >("Q_mnkl")),
  _Qmnkl(),
  _qijkl(),
  _T(coupledValue("T")),
  _T0(getParam<Real>("T0")),
  _thermal_expansion_coeff(getParam<Real>("thermal_expansion_coeff")),
  _applied_strain_vector(getParam<std::vector<Real> >("applied_strain_vector"))
//  _Electrostrictive_Euler_angles(getParam<Real>("electrostrictive_euler_angle_1"),
//                                 getParam<Real>("electrostrictive_euler_angle_2"),
//                                 getParam<Real>("electrostrictive_euler_angle_3")),

{
  //Initialize applied strain tensor from input vector
  if (_applied_strain_vector.size() == 6)
    _applied_strain_tensor.fillFromInputVector(_applied_strain_vector);
  else
    _applied_strain_tensor.zero();

  //Fill the electrostrictive coefficients
  _Qmnkl.fillFromInputVector(_Qmnkl_vector,_fill_method);

}

void
LinearFerroelectricMaterial::computeQpStrain()
{
  //strain = (grad_disp + grad_disp^T)/2
  RankTwoTensor grad_tensor(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);

  _elastic_strain[_qp] = (grad_tensor + grad_tensor.transpose()) / 2.0;
  _total_strain[_qp] = _elastic_strain[_qp];
}

void
LinearFerroelectricMaterial::computeQpStress()
{
  //Calculation and Apply stress free strain
  RankTwoTensor stress_free_strain = computeStressFreeStrain();

  // add the stress free strain on here
  // ther derivatives of elastic_strain w.r.t. c are built down in EigenstrainBaseMaterial
  _elastic_strain[_qp] += stress_free_strain;

  _total_strain[_qp] = _elastic_strain[_qp];

  // stress = C * e
  _stress[_qp] = _elasticity_tensor[_qp] * _elastic_strain[_qp];
}

RankTwoTensor
LinearFerroelectricMaterial::computeStressFreeStrain()
{
  //Apply thermal expansion
  RankTwoTensor stress_free_strain;
  stress_free_strain.addIa(-_thermal_expansion_coeff * (_T[_qp] - _T0));

  //Apply uniform applied strain
  if (_applied_strain_vector.size() == 6)
    stress_free_strain += _applied_strain_tensor;

  return stress_free_strain;
}

void
LinearFerroelectricMaterial::computeQpElectrostrictiveCoefficients()
{
// eR type: RealTensorValue
// RotationTensor eR(_Electrostrictive_Euler_angles);

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
