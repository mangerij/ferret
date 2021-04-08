//Holds temperature value on a block
#include "zT_aux.h"

registerMooseObject("FerretApp", zT_aux);

template <>
InputParameters validParams<zT_aux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates thermoelectric figure of merit");
  params.addCoupledVar("T", "Temperature variable in Kelvin");

  return params;
}


zT_aux::zT_aux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _T_var(coupled("T")),
  _T(coupledValue("T")),
  _T_grad(coupledGradient("T")),
  _ecC(getMaterialProperty<Real>("ecC")),
  _sbC(getMaterialProperty<Real>("sbC")),
  _thC(getMaterialProperty<Real>("thC"))
{
}

Real
zT_aux::computeValue()

{
    return _ecC[_qp] * _sbC[_qp] * _sbC[_qp]  / _thC[_qp] * _T[_qp];
}
