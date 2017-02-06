/************************************************************************/
/* ScreenedBC:                                                          */
/*     This BC is intended to implement a \lambda P                     */
/*     + \epsilon * _grad\phi = 0 at a boundary where                   */        
/*     \lambda is a screening param.                                    */
/************************************************************************/

#include "ScreenedBC.h"

template<>
InputParameters validParams<ScreenedBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
    params.addRequiredParam<Real>("permittivity", "permittivity");
    params.addRequiredParam<Real>("lambda", "lambda");
    params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
    params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
    params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

ScreenedBC::ScreenedBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _potential_int_grad(coupledGradient("potential_int")),
  _permittivity(getParam<Real>("permittivity")),
  _lambda(getParam<Real>("lambda")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{
}

Real
ScreenedBC::computeQpResidual()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  return _normals[_qp] * (_permittivity * _potential_int_grad[_qp] + _lambda * w);
}
