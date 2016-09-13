/**
 * @file   KappaTDiffusion.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * Note we hard-code the kappa T dependence here for PTO.
 */

#include "KappaTDiffusion.h"

template<>
InputParameters validParams<KappaTDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("temperature", "The local temperature");
  params.addParam<Real>("c0", 0.0, "the first fit coefficient");
  params.addParam<Real>("c1", 0.0, "the second fit coefficient");
  params.addParam<Real>("c2", 0.0, "the third fit coefficient");
  params.addParam<Real>("c3", 0.0, "the fourth fit coefficient");
  params.addParam<Real>("c4", 0.0, "the fifth fit coefficient");
//  params.addParam<Real>("c5", 0.0, "the sixth fit coefficient");
  return params;
}

KappaTDiffusion::KappaTDiffusion(const InputParameters & parameters)
  :Kernel(parameters),
   _temperature(coupledValue("temperature")),
   _temperature_grad(coupledGradient("temperature")),
   _c0(getParam<Real>("c0")),
   _c1(getParam<Real>("c1")),
   _c2(getParam<Real>("c2")),
   _c3(getParam<Real>("c3")),
   _c4(getParam<Real>("c4"))
{
}

Real
KappaTDiffusion::computeQpResidual()
{
  return (_c0 * _temperature[_qp] *_temperature[_qp] * std::exp(-_c1 *_temperature[_qp])+_c2*_temperature[_qp]*std::exp(- _c3 *_temperature[_qp])  + _c4 * _temperature[_qp] ) * _temperature_grad[_qp] * _grad_test[_i][_qp];
}

Real
KappaTDiffusion::computeQpJacobian()
{
  return (_c4 + _c2 *std::exp(-_c3 *_temperature[_qp]) + 2.0 * _c0 * std::exp(-_c1 *_temperature[_qp]) *_temperature[_qp] - _c2 *_c3 *std::exp(-_c3 *_temperature[_qp]) *_temperature[_qp] -  _c0* _c1 *std::exp(-_c1 *_temperature[_qp])* _temperature[_qp] *_temperature[_qp])*_phi[_j][_qp] *_temperature_grad[_qp]* _grad_test[_i][_qp] + _grad_phi[_j][_qp]* _grad_test[_i][_qp] *(_c0 * _temperature[_qp] *_temperature[_qp] * std::exp(-_c1 *_temperature[_qp])+_c2*_temperature[_qp]*std::exp(- _c3 *_temperature[_qp])  + _c4 * _temperature[_qp]);
}
