#include "DivP.h"
/// implements div of P auxkernel
template<>

InputParameters validParams<DivP>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}


DivP::DivP(const InputParameters & parameters) :
  AuxKernel(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z"))
{
}

Real
DivP::computeValue()

{
    return _polar_x_grad[_qp](0) + _polar_y_grad[_qp](1) + _polar_z_grad[_qp](2);
}
