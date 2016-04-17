/**
 * @file   CoupledEnergy.C
 * @author J. Mangeri <mangerij@anl.gov>
 *
 */

#include "CoupledEnergy.h"
#include "ComputeElectrostrictiveTensor.h"

template<>
InputParameters validParams<CoupledEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("disp_x", "The x component of the elasticity displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elasticity displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elasticity displacement");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

CoupledEnergy::CoupledEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _electrostrictive_tensor(getMaterialProperty<ElectrostrictiveTensorR4>("electrostrictive_tensor")),
  _disp_x_grad(coupledGradient("disp_x")),
  _disp_y_grad(coupledGradient("disp_y")),
  _disp_z_grad(coupledGradient("disp_z")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
CoupledEnergy::computeQpIntegral()
{
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  Real sum3 = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  sum1 += _electrostrictive_tensor[_qp].electrostrictiveProduct(0, _disp_x_grad[_qp], 0, w);
  sum1 += _electrostrictive_tensor[_qp].electrostrictiveProduct(1, _disp_y_grad[_qp], 0, w);
  sum1 += _electrostrictive_tensor[_qp].electrostrictiveProduct(2, _disp_z_grad[_qp], 0, w);

  sum2 += _electrostrictive_tensor[_qp].electrostrictiveProduct(0, _disp_x_grad[_qp], 1, w);
  sum2 += _electrostrictive_tensor[_qp].electrostrictiveProduct(1, _disp_y_grad[_qp], 1, w);
  sum2 += _electrostrictive_tensor[_qp].electrostrictiveProduct(2, _disp_z_grad[_qp], 1, w);

  sum3 += _electrostrictive_tensor[_qp].electrostrictiveProduct(0, _disp_x_grad[_qp], 2, w);
  sum3 += _electrostrictive_tensor[_qp].electrostrictiveProduct(1, _disp_y_grad[_qp], 2, w);
  sum3 += _electrostrictive_tensor[_qp].electrostrictiveProduct(2, _disp_z_grad[_qp], 2, w);

  return - 0.5 * std::pow(_len_scale, 3.0) * ( sum1 * _polar_x[_qp] + sum2 * _polar_y[_qp] + sum3 * _polar_z[_qp]);
}
