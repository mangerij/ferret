/**
 * @file   KarmanenkoDriver.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * named after Karmanenko et al J. Euro. Ceram. Soc. 27 (2007) 3109–3112
 * this term drives the temperature changes due to the field-induced entropic changes
 * NOTE: this is just a test kernel for now, as the anisotropy of _grad_potential_int
 * needs to be taken into account
 *
 * The procedure is as follows, dEstep will be related to the stepping procedure in the
 * quasi-static hysteresis loop. The only difficulty will pinning down how noise 
 * introduced is related to this kernel.
 *
 * Currently the kernel is setup for just a field along z and using the approximation of 
 * Gu et al Appl. Phys. Lett. 102, 112901, (2013).
*/

#include "KarmanenkoDriver.h"
#include<cmath>

template<>
InputParameters validParams<KarmanenkoDriver>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addCoupledVar("potential_ext", 0.0, "The external electric potential variable");
  params.addRequiredCoupledVar("temperature", "The temperature at the grid point");
  params.addParam<Real>("rho1", 1.0, "the density of the electrocaloric material");
  params.addParam<Real>("C1", 1.0, "the first coefficient of the residual contribution");
  params.addParam<Real>("C2", 1.0, "the second coefficient of the residual contribution");
  params.addParam<Real>("C3", 1.0, "the third coefficient of the residual contribution");
  params.addParam<Real>("C4", 1.0, "the fourth coefficient of the residual contribution");
  params.addParam<Real>("dEstep", 0.0, "change in the electric field as a function of time");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

KarmanenkoDriver::KarmanenkoDriver(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_ext_var(coupled("potential_ext")),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext")),
   _temperature_var(coupled("temperature")),
   _temperature(coupledValue("temperature")),
   _rho1(getParam<Real>("rho1")),
   _C1(getParam<Real>("C1")),
   _C2(getParam<Real>("C2")),
   _C3(getParam<Real>("C3")),
   _C4(getParam<Real>("C4")),
   _dEstep(getParam<Real>("dEstep")),
   _len_scale(getParam<Real>("len_scale"))
{
  std::cout<<"Implementing Karmanenko field-induced entropic change step with dEstep = "<< _dEstep <<"\n";
}

Real
KarmanenkoDriver::computeQpResidual()

///0.000116379 Ez + 0.00117004 Ez^2 - 3.27924*10^-6 Ez T + 4.68114*10^-6 Ez^2 T
{
  return _test[_i][_qp] * std::pow(_len_scale, 2.0) * _rho1 * _temperature[_qp] * (- _C1 * _potential_int_grad[_qp](2) + _C2 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) - _C3 * _potential_int_grad[_qp](2) * _temperature[_qp] + _C4 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) * _temperature[_qp]) * _dEstep;
}

Real
KarmanenkoDriver::computeQpJacobian()
{
  return _test[_i][_qp] * std::pow(_len_scale, 2.0) * _rho1 * _phi[_j][_qp] * (- _C1 * _potential_int_grad[_qp](2) + _C2 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) - _C3 * _potential_int_grad[_qp](2) * _temperature[_qp] + _C4 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) * _temperature[_qp]) * _dEstep;
}

Real
KarmanenkoDriver::computeQpOffDiagJacobian(unsigned int jvar)
{
    if( jvar == _potential_int_var )
      return  _test[_i][_qp] * std::pow(_len_scale, 2.0) * _rho1 * _temperature[_qp] * (- _C1 * _grad_phi[_j][_qp](2) + _C2 * _grad_phi[_j][_qp](2) * _potential_int_grad[_qp](2) - _C3 * _grad_phi[_j][_qp](2) * _temperature[_qp] + _C4 * _grad_phi[_j][_qp](2) * _potential_int_grad[_qp](2) * _temperature[_qp]) * _dEstep;;
    else if( jvar == _potential_ext_var)
      return  0.0;
    else
    {
      return 0.0;
    }
}
