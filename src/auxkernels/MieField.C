#include "MieField.h"
#include <complex>
#include <cmath>

template<>

InputParameters validParams<MieField>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("a", "radius of the particle");
  params.addRequiredParam<Real>("omega", "angular frequency of the plane wave");
  params.addRequiredParam<Real>("c", "the speed of light");
  params.addRequiredParam<Real>("epsilonII", "the dielectric constant of the sphere");
  params.addRequiredParam<Real>("sigmaII", "the conductivity of the sphere");
  params.addRequiredParam<Real>("L", "wavelength of the plane wave");
  params.addRequiredParam<Real>("order", "order of the expansions");
  params.addRequiredParam<Real>("component", "component at which to compute this field ('0 for x, 1 for y, 2 for z')");
  return params;
}


MieField::MieField(const InputParameters & parameters) :
  AuxKernel(parameters),
  _a(getParam<Real>("a")),
  _omega(getParam<Real>("omega")),
  _c(getParam<Real>("c")),
  _epsilonII(getParam<Real>("epsilonII")),
  _sigmaII(getParam<Real>("sigmaII")),
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

  Complex k1 = Complex(0, 1) * (_omega / _c) * (_epsilonII + Complex(0, 1) * 4.0 * libMesh::pi * _sigmaII / _omega);
  Complex k2 = Complex(0, 1) * (_omega / _c);

  Complex ksq = - k1 * k2;
  
  Real r =  std::pow(x, 2.0) + std::pow(y, 2.0) + std::pow(z, 2.0);
  Real th = std::acos(z / r);
  Real phi = std::atan(y / x);
  Real q = 2 * libMesh::pi * _a / _L;

  Complex Ersj = 0.0;
  Complex Ethsj = 0.0;
  Complex Ephisj = 0.0;
 
  Complex Ex = 0.0;
  Complex Ey = 0.0;
  Complex Ez = 0.0;

  for(unsigned int j = 1; j < _order; ++j)
    {
       Complex eBnumerator = 1.0;
       Complex mBnumerator = 1.0;
       Complex eBdenominator = 1.0;
       Complex mBdenominator = 1.0;
       Complex eBj = ( (2 * j + 1) / (j + 1) ); //* eBnumerator / eBdenominator;
       Complex mBj = ( (2 * j + 1) / (j + 1) ); //* mBnumerator / mBdenominator;
       eBj = eBj * Complex(0,1); //this is likely what one needs to do...
       Ersj += (std::pow(_L, 2.0) * std::cos(phi) / (4 * std::pow(libMesh::pi, 2.0) * r * r)) * j * (j + 1) * eBj;
       Ethsj += (std::pow(_L, 2.0) * std::cos(phi)/ (4 * std::pow(libMesh::pi, 2.0) * r)) * eBj * std::sin(th);
       Ephisj += (std::pow(_L, 2.0) * std::sin(phi) / (4 * std::pow(libMesh::pi, 2.0) * r)) * eBj / std::sin(th);
    }

  if (_component == 0)
    {
      Ex += Ethsj * std::cos(th) * std::cos(phi) +  Ersj * std::cos(phi) * std::sin(th) - Ephisj * std::sin(phi);
      return Ex.real();
    }
  else if (_component == 1)
    {
      Ey += Ersj * std::sin(th) * std::sin(phi) +  std::cos(phi) *(Ephisj + Ethsj * std::sin(phi) );
      return Ey.real();
    }
  else if (_component == 2)
    {
      Ez += Ersj * std::cos(th) - Ethsj * std::sin(phi);
      return Ez.real();
    }
  else
    {
      return 0.0;
    }
}


