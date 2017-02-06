/**
 * @file   DepolarizationEnergy.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief
 *
 *
 */

#include "DepolarizationEnergy.h"

template<>
InputParameters validParams<DepolarizationEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("avePz", 1.0, "the average polarization due to the depol field term");
  params.addParam<Real>("lambda", 1.0, "the screening length term");
  params.addParam<Real>("permitivitty", 1.0, "the permitivitty term");
  return params;
}

DepolarizationEnergy::DepolarizationEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _avePz(getParam<Real>("avePz")),
   _lambda(getParam<Real>("lambda")),
   _permitivitty(getParam<Real>("permitivitty"))
{
}

Real
DepolarizationEnergy::computeQpIntegral()
{
  return 0.5 * _lambda * (1.0 / _permitivitty) * _polar_z[_qp] * _avePz * std::pow(_len_scale, 3.0);
}
