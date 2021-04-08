#include "ResidualT.h"
#include "Material.h"

registerMooseObject("FerretApp", ResidualT);

template <>
InputParameters
validParams<ResidualT>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to modified ohm's law");
  params.addRequiredCoupledVar("potential_E_int", "electrical potential");
  params.addRequiredCoupledVar("T", "temperature");
  return params;
}

ResidualT::ResidualT(const InputParameters & parameters)
  : Kernel(parameters),
    _potential_E_int_var(coupled("potential_E_int")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _T_var(coupled("T")),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _ecC(getMaterialProperty<Real>("ecC")),
    _sbC(getMaterialProperty<Real>("sbC")),
    _thC(getMaterialProperty<Real>("thC"))
{
}

Real
ResidualT::computeQpResidual()
{
  return ((-_grad_test[_i][_qp](0)) * (-_thC[_qp] * _T_grad[_qp](0)) +
          (-_grad_test[_i][_qp](0)) *
              (-_ecC[_qp] * _sbC[_qp] * _sbC[_qp] * _T[_qp] * _T_grad[_qp](0)) +
          (-_grad_test[_i][_qp](0)) *
              (-_ecC[_qp] * _sbC[_qp] * _T[_qp] * _potential_E_int_grad[_qp](0)) +
          _test[_i][_qp] *
              (-_ecC[_qp] * _potential_E_int_grad[_qp](0) * _potential_E_int_grad[_qp](0)) +
          _test[_i][_qp] *
              (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](0) * _potential_E_int_grad[_qp](0))) +
         ((-_grad_test[_i][_qp](1)) * (-_thC[_qp] * _T_grad[_qp](1)) +
          (-_grad_test[_i][_qp](1)) *
              (-_ecC[_qp] * _sbC[_qp] * _sbC[_qp] * _T[_qp] * _T_grad[_qp](1)) +
          (-_grad_test[_i][_qp](1)) *
              (-_ecC[_qp] * _sbC[_qp] * _T[_qp] * _potential_E_int_grad[_qp](1)) +
          _test[_i][_qp] *
              (-_ecC[_qp] * _potential_E_int_grad[_qp](1) * _potential_E_int_grad[_qp](1)) +
          _test[_i][_qp] *
              (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](1) * _potential_E_int_grad[_qp](1))) +
         ((-_grad_test[_i][_qp](2)) * (-_thC[_qp] * _T_grad[_qp](2)) +
          (-_grad_test[_i][_qp](2)) *
              (-_ecC[_qp] * _sbC[_qp] * _sbC[_qp] * _T[_qp] * _T_grad[_qp](2)) +
          (-_grad_test[_i][_qp](2)) *
              (-_ecC[_qp] * _sbC[_qp] * _T[_qp] * _potential_E_int_grad[_qp](2)) +
          _test[_i][_qp] *
              (-_ecC[_qp] * _potential_E_int_grad[_qp](2) * _potential_E_int_grad[_qp](2)) +
          _test[_i][_qp] *
              (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](2) * _potential_E_int_grad[_qp](2)));
}

Real
ResidualT::computeQpJacobian()
{
  return ((-_grad_test[_i][_qp](0) *
               (-_ecC[_qp] * _sbC[_qp] * _phi[_j][_qp] * _potential_E_int_grad[_qp](0)) +
           (-_grad_test[_i][_qp](0)) * (-_thC[_qp] * _grad_phi[_j][_qp](0)) +
           (-_grad_test[_i][_qp](0)) *
               (-_sbC[_qp] * _sbC[_qp] * _ecC[_qp] *
                (_phi[_j][_qp] * _T_grad[_qp](0) + _T[_qp] * _grad_phi[_j][_qp](0))) +
           _test[_i][_qp] *
               (-_sbC[_qp] * _ecC[_qp] * _grad_phi[_j][_qp](0) * _potential_E_int_grad[_qp](0))) +
          ((-_grad_test[_i][_qp](1)) *
               (-_ecC[_qp] * _sbC[_qp] * _phi[_j][_qp] * _potential_E_int_grad[_qp](1)) +
           (-_grad_test[_i][_qp](1)) * (-_thC[_qp] * _grad_phi[_j][_qp](1)) +
           (-_grad_test[_i][_qp](1)) *
               (-_sbC[_qp] * _sbC[_qp] * _ecC[_qp] *
                (_phi[_j][_qp] * _T_grad[_qp](1) + _T[_qp] * _grad_phi[_j][_qp](1))) +
           _test[_i][_qp] *
               (-_sbC[_qp] * _ecC[_qp] * _grad_phi[_j][_qp](1) * _potential_E_int_grad[_qp](1))) +
          ((-_grad_test[_i][_qp](2)) *
               (-_ecC[_qp] * _sbC[_qp] * _phi[_j][_qp] * _potential_E_int_grad[_qp](2)) +
           (-_grad_test[_i][_qp](2)) * (-_thC[_qp] * _grad_phi[_j][_qp](2)) +
           (-_grad_test[_i][_qp](2)) *
               (-_sbC[_qp] * _sbC[_qp] * _ecC[_qp] *
                (_phi[_j][_qp] * _T_grad[_qp](2) + _T[_qp] * _grad_phi[_j][_qp](2))) +
           _test[_i][_qp] *
               (-_sbC[_qp] * _ecC[_qp] * _grad_phi[_j][_qp](2) * _potential_E_int_grad[_qp](2))));
}

Real
ResidualT::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_E_int_var)
  {
    return (
        (-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _sbC[_qp] * _T[_qp] * _grad_phi[_j][_qp](0)) -
        _test[_i][_qp] * (2.0 * _ecC[_qp] * _potential_E_int_grad[_qp](0) * _grad_phi[_j][_qp](0) +
                          _ecC[_qp] * _sbC[_qp] * _T_grad[_qp](0) * _grad_phi[_j][_qp](0)) +
        (-_grad_test[_i][_qp](1)) * (-_ecC[_qp] * _sbC[_qp] * _T[_qp] * _grad_phi[_j][_qp](1)) -
        _test[_i][_qp] * (2.0 * _ecC[_qp] * _potential_E_int_grad[_qp](1) * _grad_phi[_j][_qp](1) +
                          _ecC[_qp] * _sbC[_qp] * _T_grad[_qp](1) * _grad_phi[_j][_qp](1)) +
        (-_grad_test[_i][_qp](2)) * (-_ecC[_qp] * _sbC[_qp] * _T[_qp] * _grad_phi[_j][_qp](2)) -
        _test[_i][_qp] * (2.0 * _ecC[_qp] * _potential_E_int_grad[_qp](2) * _grad_phi[_j][_qp](2) +
                          _ecC[_qp] * _sbC[_qp] * _T_grad[_qp](2) * _grad_phi[_j][_qp](2)));
  }
  else
    return 0.0;
}
