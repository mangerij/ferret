#include "SemiconductingChargeCarriersAux.h"

template<>

InputParameters validParams<SemiconductingChargeCarriersAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("charge_type", "An integer corresponding to the charge type to calculate (0 for nm, 1 for pp, 2 for NAm, 3 for total rho)");
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential");
  params.addParam<Real>("q", 1.0, "the variable q corresponding to the value of the charge");
  params.addParam<Real>("kT", 1.0, "kT");
  params.addParam<Real>("NA", 1.0, "NA");
  params.addParam<Real>("NC", 1.0, "NC");
  params.addParam<Real>("NV", 1.0, "NV");
  params.addParam<Real>("EA", 1.0, "EA");
  params.addParam<Real>("EC", 1.0, "EC");
  params.addParam<Real>("EV", 1.0, "EV");
  params.addParam<Real>("EF", 1.0, "EF");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}


SemiconductingChargeCarriersAux::SemiconductingChargeCarriersAux(const InputParameters & parameters) :
  AuxKernel(parameters ),
   _charge_type(getParam<unsigned int>("charge_type")),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _q(getParam<Real>("q")),
   _kT(getParam<Real>("kT")),
   _NA(getParam<Real>("NA")),
   _NC(getParam<Real>("NC")),
   _NV(getParam<Real>("NV")),
   _EA(getParam<Real>("EA")),
   _EC(getParam<Real>("EC")),
   _EV(getParam<Real>("EV")),
   _EF(getParam<Real>("EF"))
{
}

Real
SemiconductingChargeCarriersAux::computeValue()

{
  if (_charge_type == 0)
  {
    Real nm = _NC * std::pow((std::exp(_q * (_EC - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0);
    return nm;
  }
  else if (_charge_type == 1)
  {
    Real pp = _NV * (1 - std::pow((std::exp(_q * (_EV - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0));
    return pp;
  }
  else if (_charge_type == 2)
  {
    Real NAm = _NA * std::pow((std::exp(_q * (_EA - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0);
    return NAm;
  }
  else if (_charge_type == 3)
  {
    Real nm = _NC * std::pow((std::exp(_q * (_EC - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0);
    Real pp = _NV * (1 - std::pow((std::exp(_q * (_EV - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0));
    Real NAm = _NA * std::pow((std::exp(_q * (_EA - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0);
    Real rho = _q * ( - nm + pp - NAm);
    return rho;
  }
else 
  return 0.0;
}
