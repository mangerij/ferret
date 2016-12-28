/**
 * @file   CubicDielectricTensor.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief calculate the components of the anisotropic dielectric tensor
 *        assuming cubic symmetry of the parent phase
 *
 */
#include "CubicDielectricTensor.h"
#include "ComputeElectrostrictiveTensor.h"

template<>
InputParameters validParams<CubicDielectricTensor>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0,  "The z component of the elastic displacement");
  params.addRequiredParam<Real>("alpha1", "alpha1 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "alpha11 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "alpha12 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "alpha111 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "alpha112 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "alpha123 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("first_deriv", "direction of first derivative");
  params.addRequiredParam<Real>("second_deriv", "direction of second derivative");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

CubicDielectricTensor::CubicDielectricTensor(const InputParameters & parameters) :
  AuxKernel(parameters),
  _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _disp_x_grad(coupledGradient("disp_x")),
  _disp_y_grad(coupledGradient("disp_y")),
  _disp_z_grad(coupledGradient("disp_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha12(getParam<Real>("alpha12")),
  _alpha111(getParam<Real>("alpha111")),
  _alpha112(getParam<Real>("alpha112")),
  _alpha123(getParam<Real>("alpha123")),
  _first_deriv(getParam<Real>("first_deriv")),
  _second_deriv(getParam<Real>("second_deriv")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
CubicDielectricTensor::computeValue()
{
  //normal components:
  if (_first_deriv == 0 && _second_deriv == 0)
  {
    return 0.0;
    // return 2 * (
    //    _alpha1 + 15.0 * std::pow(_polar_x[_qp], 4.0) * _alpha111 + _alpha112 * std::pow(_polar_y[_qp], 4.0)
    //    + _alpha112 * std::pow(_polar_z[_qp], 4.0) + 6.0 * _polar_x[_qp] *_polar_x[_qp] * (_alpha11
    //   + _alpha112 * (_polar_y[_qp] *_polar_y[_qp] + _polar_z[_qp] *_polar_z[_qp])
    // ) + _polar_z[_qp] * _polar_z[_qp] * _alpha12 + _polar_y[_qp] * _polar_y[_qp] * (_alpha12 + _polar_z[_qp] * _polar_z[_qp] * _alpha123)
    // + _electrostrictive_tensor[_qp](0, 0, 0, 0) * _disp_x_grad[_qp](0) + 0.5 * _electrostrictive_tensor[_qp](0, 0, 1, 1) * (_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0))
    // + 0.5 * _electrostrictive_tensor[_qp](0, 0, 2, 2) * (_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)));
  }
  else if (_first_deriv == 1 && _second_deriv == 1)
  {
    return 0.0;
    // return 2*(_alpha1 + 15*std::power(_polar_y[_qp], 4.0)*_alpha111 + std::power(_polar_x[_qp],4)*_alpha112 + 6*std::power(_polar_y[_qp], 2.0)*(_alpha11
    //   + (std::power(_polar_x[_qp], 2.0) + std::power(Pz,2))*_alpha112) + std::power(Px,2)*_alpha12
    //   + std::power(Pz,2)*(std::power(Pz,2)*_alpha112 + _alpha12 + std::power(Px,2)*_alpha123) +
    //   _electrostrictive_tensor[_qp](0,0,1,1)*epsilon(0,1) + _electrostrictive_tensor[_qp](1,1,1,1)*epsilon(1,1) + _electrostrictive_tensor[_qp](2,1)*epsilon(2,2));
  }
  else if (_first_deriv == 2 && _second_deriv == 2)
  {
    return 0.0;
  }
  //shears:
  else if (_first_deriv == 0 && _second_deriv == 1)
  {
    return 0.0;
  }
  else if (_first_deriv == 0 && _second_deriv == 2)
  {
    return 0.0;
  }
  else if (_first_deriv == 1 && _second_deriv == 2)
  {
    return 0.0;
  }
  else
    return 0.0;
}
