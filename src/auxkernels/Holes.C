#include "Holes.h"
registerMooseObject("FerretApp", Holes);

InputParameters Holes::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<Real>("Ev", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nv", "Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredCoupledVar("potential_E_int","E");
  return params;
}


Holes::Holes(const InputParameters & parameters) :
  AuxKernel(parameters),
  _Ev(getParam<Real>("Ev")),
  _Nv(getParam<Real>("Nv")),
  _T(getParam<Real>("T")),
  _Kb(getParam<Real>("Kb")),
  _q(getParam<Real>("q")),
  _potential_E_int(coupledValue("potential_E_int"))
{
}

Real
Holes::computeValue()

{
    return (_Nv * std::exp((( _Ev - (_q * _potential_E_int[_qp])) / (_Kb * _T))));
}
