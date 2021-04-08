#include "ResidualV.h"
// #include "HeatConduction.h"
#include "Material.h"

registerMooseObject("FerretApp", ResidualV);

template <>
InputParameters
validParams<ResidualV>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to modified ohm's law");
  params.addRequiredCoupledVar("potential_E_int", "electrical potential");
  params.addRequiredCoupledVar("T", "temperature");
  return params;
}

ResidualV::ResidualV(const InputParameters & parameters)
  : Kernel(parameters),
    _potential_E_int_var(coupled("potential_E_int")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _T_var(coupled("T")),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _ecC(getMaterialProperty<Real>("ecC")),
    _sbC(getMaterialProperty<Real>("sbC"))
{
}

Real
ResidualV::computeQpResidual()

{
  return (((-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](0)) +
           (-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _potential_E_int_grad[_qp](0))) +
          (-_grad_test[_i][_qp](1) * (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](1)) +
           (-_grad_test[_i][_qp](1)) * (-_ecC[_qp] * _potential_E_int_grad[_qp](1))) +
          (-_grad_test[_i][_qp](2) * (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](2)) +
           (-_grad_test[_i][_qp](2)) * (-_ecC[_qp] * _potential_E_int_grad[_qp](2))));
}

Real
ResidualV::computeQpJacobian()

{
  return ((-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _grad_phi[_j][_qp](0)) +
          (-_grad_test[_i][_qp](1)) * (-_ecC[_qp] * _grad_phi[_j][_qp](1)) +
          (-_grad_test[_i][_qp](2)) * (-_ecC[_qp] * _grad_phi[_j][_qp](2)));
}

Real
ResidualV::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _T_var)
  {
    return ((-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _sbC[_qp] * _grad_phi[_j][_qp](0)) +
            (-_grad_test[_i][_qp](1)) * (-_ecC[_qp] * _sbC[_qp] * _grad_phi[_j][_qp](1)) +
            (-_grad_test[_i][_qp](2)) * (-_ecC[_qp] * _sbC[_qp] * _grad_phi[_j][_qp](2)));
  }
  else
  {
    return 0.0;
  }
}
