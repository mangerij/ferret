/**
 * @file   CoupledEnergy.C
 * @author J. Mangeri <mangerij@anl.gov>
 *
 */

#include "CoupledEnergy.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ComputeEigenstrain.h"

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
  params.addParam<Real>("artificial", 1.0, "artificial increase coupling");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

CoupledEnergy::CoupledEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
  _stress_free_strain(getMaterialProperty<RankTwoTensor>("stress_free_strain")),
  _disp_x_grad(coupledGradient("disp_x")),
  _disp_y_grad(coupledGradient("disp_y")),
  _disp_z_grad(coupledGradient("disp_z")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _artificial(getParam<Real>("artificial")),
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
  RealVectorValue v0(_stress_free_strain[_qp](0,0), _stress_free_strain[_qp](0,1), _stress_free_strain[_qp](0,2));
  RealVectorValue v1(_stress_free_strain[_qp](1,0), _stress_free_strain[_qp](1,1), _stress_free_strain[_qp](1,2));
  RealVectorValue v2(_stress_free_strain[_qp](2,0), _stress_free_strain[_qp](2,1), _stress_free_strain[_qp](2,2));

  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, 0, w);
  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, 0, w);
  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, 0, w);

  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, 1, w);
  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, 1, w);
  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, 1, w);

  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, 2, w);
  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, 2, w);
  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, 2, w);

  return -_artificial * std::pow(_len_scale, 3.0) * ( sum1 * _polar_x[_qp] + sum2 * _polar_y[_qp] + sum3 * _polar_z[_qp]);
}
