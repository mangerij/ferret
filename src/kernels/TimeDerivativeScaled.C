
/* Modified(Scaled) for FERRET, added time_scale  */

#include "TimeDerivativeScaled.h"

// MOOSE includes
#include "MooseVariable.h"

template<>
InputParameters validParams<TimeDerivativeScaled>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  params.addParam<Real>("time_scale",1.0,"the time_scale of the unit");
  return params;
}

TimeDerivativeScaled::TimeDerivativeScaled(const InputParameters & parameters) :
    TimeKernel(parameters),
    _lumping(getParam<bool>("lumping")),
    _time_scale(getParam<Real>("time_scale"))
{
}

Real
TimeDerivativeScaled::computeQpResidual()
{
  return _time_scale * _test[_i][_qp] * _u_dot[_qp];
}

Real
TimeDerivativeScaled::computeQpJacobian()
{
  return _time_scale * _test[_i][_qp]*_phi[_j][_qp] * _du_dot_du[_qp];
}

void
TimeDerivativeScaled::computeJacobian()
{
  if (_lumping)
  {
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
          ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
        }
  }
  else
    TimeKernel::computeJacobian();
}
