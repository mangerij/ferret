#include "HoleGen.h"

registerMooseObject("FerretApp", HoleGen);

template<>
InputParameters validParams<HoleGen>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("Ev", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nv", "Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("mup","hole mobility");
  return params;
}

HoleGen::HoleGen(const InputParameters & parameters)
  :Kernel(parameters),
  _Ev(getParam<Real>("Ev")),
  _Nv(getParam<Real>("Nv")),
  _T(getParam<Real>("T")),
  _Kb(getParam<Real>("Kb")),
  _q(getParam<Real>("q")),
  _mup(getParam<Real>("mup"))
{
}

Real
HoleGen::computeQpResidual()
{
  Real Rhgen = 0.0;
  Rhgen += (_mup * _Kb * _T / _q) *

  ((_grad_test[_i][_qp](0) * _Nv * std::exp(( _Ev - (_q * _grad_u[_qp](0))) / (_Kb * _T))) +

  (_grad_test[_i][_qp](1) * _Nv * std::exp(( _Ev - (_q * _grad_u[_qp](1))) / (_Kb * _T))) +

  (_grad_test[_i][_qp](2) * _Nv * std::exp(( _Ev - (_q * _grad_u[_qp](2))) / (_Kb * _T))));

  return Rhgen;
}

Real
HoleGen::computeQpJacobian()
{
   return (_mup * _Kb * _T / _q) *

   (((_grad_test[_i][_qp](0) * -1 * _Nv * _q * _grad_phi[_j][_qp](0) / (_Kb * _T)) * std::exp(( _Ev - (_q * _grad_u[_qp](0))) / (_Kb * _T))) +

   (( _grad_test[_i][_qp](1) * -1 * _Nv * _q * _grad_phi[_j][_qp](1) / (_Kb * _T)) * std::exp(( _Ev - (_q * _grad_u[_qp](1))) / (_Kb * _T))) +

   (( _grad_test[_i][_qp](2) * -1 * _Nv * _q * _grad_phi[_j][_qp](2) / (_Kb * _T)) * std::exp(( _Ev - (_q * _grad_u[_qp](2))) / (_Kb * _T))));
}
