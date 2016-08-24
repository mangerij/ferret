#include "MieField.h"

template<>

InputParameters validParams<MieField>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("a", "radius of the particle");
  params.addRequiredParam<Real>("L", "wavelength of the plane wave");
  params.addRequiredParam<Real>("order", "order of the expansions");
  params.addRequiredParam<Real>("component", "component at which to compute this field ('0 for x, 1 for y, 2 for z')");
  return params;
}


MieField::MieField(const InputParameters & parameters) :
  AuxKernel(parameters),
  _a(getParam<Real>("a")),
  _L(getParam<Real>("L")),
  _order(getParam<Real>("order")),
  _component(getParam<Real>("component"))
{
}

Real
MieField::computeValue()
{
  Real x = _q_point[_qp](0);
  Real y = _q_point[_qp](1);
  Real z = _q_point[_qp](2);

  Real r =  std::pow(x, 2.0) + std::pow(y, 2.0) + std::pow(z, 2.0);
  Real th = std::acos(z / r);
  Real phi = std::atan(y / x);
  Real q = 2 * 3.14159 * _a / _L;

  Real Ersj = 0.0;
  Real Ethsj = 0.0;
  Real Ephisj = 0.0;
 
 // for (j = 1 to _order);
 //   {
 //     Ersj += (std::pow(_L, 2.0) * std::cos(phi) / (4 * std::pow(3.14159, 2.0) * r * r)) * j * (j + 1);
 //     Ethsj += (std::pow(_L, 2.0) * std::cos(phi)/ (4 * std::pow(3.14159, 2.0) * r));
 //     Ephisj += (std::pow(_L, 2.0) * std::sin(phi) / (4 * std::pow(3.14159, 2.0) * r));
 //   }

 // Real Ex = Ethsj * std::cos(th) * std::cos(phi) +  Ersj * std::cos(phi) * std::sin(th) - Ephisj * std::sin(phi);
 // Real Ey = Ersj * std::sin(th) * std::sin(phi) +  std::cos(phi) *(Ephisj + Ethsj * std::sin(phi) );
 // Real Ez = Ersj * std::cos(th) - Ethsj * std::sin(phi);
 return 0.0;
}


