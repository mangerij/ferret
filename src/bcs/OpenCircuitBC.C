/************************************************************************/
/* OpenCircuit BC:                                                      */
/*     This BC is intended to implement a P + \epsilon * _grad\phi = 0  */
/*     at a boundary, currently assumed to be along z.                  */
/************************************************************************/


#include "OpenCircuitBC.h"

template<>
InputParameters validParams<OpenCircuitBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
    params.addRequiredParam<Real>("permittivity", "permittivity");
    params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
    params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
    params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

OpenCircuitBC::OpenCircuitBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _potential_int_grad(coupledGradient("potential_int")),
  _permittivity(getParam<Real>("permittivity")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{
}

Real
OpenCircuitBC::computeQpResidual()
{
  return _permittivity * _potential_int_grad[_qp](2) + _polar_z[_qp];
}
