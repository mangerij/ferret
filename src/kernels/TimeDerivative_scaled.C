/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
/* Modified(_scaled) for FERRET, added time_scale for ferroic dynamics   */

#include "TimeDerivative_scaled.h"

template<>
InputParameters validParams<TimeDerivative_scaled>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  params.addParam<Real>("time_scale",1.0,"the time_scale of the unit");
  return params;
}

TimeDerivative_scaled::TimeDerivative_scaled(const std::string & name, InputParameters parameters) :
    TimeKernel(name, parameters),
    _lumping(getParam<bool>("lumping")),
    _time_scale(getParam<Real>("time_scale"))
{
}

Real
TimeDerivative_scaled::computeQpResidual()
{
  return _time_scale*_test[_i][_qp]*_u_dot[_qp];
}

Real
TimeDerivative_scaled::computeQpJacobian()
{
  return _time_scale*_test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp];
}

void
TimeDerivative_scaled::computeJacobian()
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
