/**
 * @file   AnisotropicEnergy.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief
 *   adds anisotropy energy postprocessor for skyrmion problem
 *
 */


#include "AnisotropicEnergy.h"

template<>
InputParameters validParams<AnisotropicEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("K", 1.0, "the anisotropy energy");
  return params;
}

AnisotropicEnergy::AnisotropicEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _K(getParam<Real>("K"))
{
}

Real
AnisotropicEnergy::computeQpIntegral()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  return _K * (w(0) * w(0) + w(1) * w(1)) * std::pow(_len_scale, 3.0);
}
