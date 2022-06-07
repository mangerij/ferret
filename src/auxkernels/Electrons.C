#include "Electrons.h"
registerMooseObject("FerretApp", Electrons);

InputParameters Electrons::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<Real>("Ec", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nc", "Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredCoupledVar("potential_E_int","E");
  return params;
}


Electrons::Electrons(const InputParameters & parameters) :
  AuxKernel(parameters),
  _Ec(getParam<Real>("Ec")),
  _Nc(getParam<Real>("Nc")),
  _T(getParam<Real>("T")),
  _Kb(getParam<Real>("Kb")),
  _q(getParam<Real>("q")),
  _potential_E_int(coupledValue("potential_E_int"))
{
}

Real
Electrons::computeValue()

{
    return (_Nc * std::exp((((_q * _potential_E_int[_qp])-_Ec) / (_Kb * _T))));
}
