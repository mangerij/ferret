/************************************************************************/
/* OpenCircuit BC:                                                      */
/*     This BC is intended to implement a P + \epsilon * _grad\phi = 0  */
/*     at a boundary.                                                   */
/************************************************************************/

#include "OpenCircuitBC.h"

template<>
InputParameters validParams<OpenCircuitBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
    params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
    params.addRequiredParam<Real>("permittivity", "permittivity");
    params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
    params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
    params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

OpenCircuitBC::OpenCircuitBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _component(getParam<unsigned int>("component")),
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
  const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1) ? _polar_y: _polar_z;

  return _normals[_qp](_component) * (_permittivity * _potential_int_grad[_qp](_component) + _polar_i[_qp]);
}
