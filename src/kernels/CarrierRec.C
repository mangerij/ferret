#include "CarrierRec.h"

registerMooseObject("FerretApp", CarrierRec);

template<>
InputParameters validParams<CarrierRec>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("Ev",  "Property name of the Valence band energy (J)");
  params.addRequiredParam<Real>("Ec", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nv","Effective DOS of the valence band(T=298)");
  params.addRequiredParam<Real>("Nc","Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("b", "lifetime of carriers");
  return params;
}

CarrierRec::CarrierRec(const InputParameters & parameters)
  :Kernel(parameters),
  _Ev(getParam<Real>("Ev")),
  _Ec(getParam<Real>("Ec")),
  _Nv(getParam<Real>("Nv")),
  _Nc(getParam<Real>("Nc")),
   _T(getParam<Real>("T")),
   _Kb(getParam<Real>("Kb")),
   _q(getParam<Real>("q")),
   _b(getParam<Real>("b"))
{
}

Real
CarrierRec::computeQpResidual()
{
  Real Rrec = 0.0;
  Rrec += -1 * _b * _test[_i][_qp] *

  (_Nc * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T))) * (_Nv * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T)));

  return Rrec;
}

Real
CarrierRec::computeQpJacobian()
{
   return -1 * _b * _test[_i][_qp] *

   ((((_Nc * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T)))*(_Nv * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T)))) +

   (((-1 * _Nv * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T)))*(_Nc * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T)))));
}
