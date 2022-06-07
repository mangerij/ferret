#include "CarrierCon.h"

registerMooseObject("FerretApp", CarrierCon);

template<>
InputParameters validParams<CarrierCon>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("Ev", "Valence band energy");
  params.addRequiredParam<Real>("Ec", "Conduction band energy");
  params.addRequiredParam<Real>("Nv","Effective DOS of the valence band(T=298)");
  params.addRequiredParam<Real>("Nc","Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("Nd", "Cation Donor");
  params.addRequiredParam<Real>("Na", "Anion Donor");
  return params;
}

CarrierCon::CarrierCon(const InputParameters & parameters)
  :Kernel(parameters),
   _Ev(getParam<Real>("Ev")),
   _Ec(getParam<Real>("Ec")),
   _Nv(getParam<Real>("Nv")),
   _Nc(getParam<Real>("Nc")),
   _T(getParam<Real>("T")),
   _Kb(getParam<Real>("Kb")),
   _q(getParam<Real>("q")),
   _Na(getParam<Real>("Na")),
   _Nd(getParam<Real>("Nd"))

{
}

Real
CarrierCon::computeQpResidual()
{
  Real Rcar = 0.0;
  Rcar += _test[_i][_qp] * _q *

  (((_Nc * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T))) -

  (_Nv * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T)))) +

  _Na - _Nd);

  return Rcar;
}

Real
CarrierCon::computeQpJacobian()
{
   return _test[_i][_qp] * _q *

   ((_Nc * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T)) -

   (-_Nv * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp(( _Ev - ( _q * _u[_qp])) / (_Kb * _T))) ;
}
