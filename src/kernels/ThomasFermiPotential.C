/**
 * @file   ThomasFermiPotential.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#include "ThomasFermiPotential.h"

class ThomasFermiPotential;

template<>
InputParameters validParams<ThomasFermiPotential>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential to be screened");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  params.addParam<Real>("TFconstant", 0.0, "the Thomas-Fermi constant which proportional to rho divided by the Fermi level");
  return params;
}

ThomasFermiPotential::ThomasFermiPotential(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _len_scale(getParam<Real>("len_scale")),
   _TFconstant(getParam<Real>("TFconstant"))
{
}

Real
ThomasFermiPotential::computeQpResidual()
{
  return _TFconstant * _potential_int[_qp] * _test[_i][_qp];
}
Real
ThomasFermiPotential::computeQpJacobian()
{
  return _TFconstant * _potential_int[_qp] * _phi[_j][_qp];
}

