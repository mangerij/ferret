/**
 * @file   AnisotropyEnergy.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief  Implement the kernel for polar variables corresponding to ferroelectic anisotropy energy after
 *         the variational derivative of the polar dependent terms have been taken.
 */

#include "AnisotropyEnergy.h"

class AnisotropyEnergy;

template<>
InputParameters validParams<AnisotropyEnergy>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("K", 1.0, "the anisotropy energy");
  return params;
}

AnisotropyEnergy::AnisotropyEnergy(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
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
AnisotropyEnergy::computeQpResidual()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  return 2.0 * _K * w(_component) * _test[_i][_qp];
}

Real
AnisotropyEnergy::computeQpJacobian()
{
  return 2.0 * _K * _phi[_j][_qp] * _test[_i][_qp];
}
