/**
 * @file   SemiconductorChargeCarriers.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Here this file will couple Boltzmann carrier distributions
 * to the electrostatic fields in a semiconducting block.
 */

#include "SemiconductorChargeCarriers.h"

class SemiconductorChargeCarriers;

template<>
InputParameters validParams<SemiconductorChargeCarriers>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

SemiconductorChargeCarriers::SemiconductorChargeCarriers(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
SemiconductorChargeCarriers::computeQpResidual()
{
  return 0.0;
}
Real
SemiconductorChargeCarriers::computeQpJacobian()
{
  return 0.0;
}

Real
SemiconductorChargeCarriers::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
