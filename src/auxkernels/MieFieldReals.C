#include "MieFieldReals.h"
#include <complex>
#include <cmath>

#ifdef FERRET_HAVE_BOOST_MATH_SPECIAL_FUNCTIONS
# include <boost/math/special_functions/legendre.hpp>
#endif

template<>

InputParameters validParams<MieFieldReals>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("a", "radius of the particle");
  params.addRequiredParam<Real>("omega", "angular frequency of the plane wave");
  params.addRequiredParam<Real>("c", "the speed of light");
  params.addRequiredParam<Real>("epsilonI", "the dielectric constant of the medium");
  params.addRequiredParam<Real>("sigmaI", "the conductivity of the medium");
  params.addRequiredParam<Real>("epsilonII", "the dielectric constant of the sphere");
  params.addRequiredParam<Real>("sigmaII", "the conductivity of the sphere");
  params.addRequiredParam<Real>("L", "wavelength of the plane wave");
  params.addRequiredParam<Real>("order", "order of the expansions");
  params.addRequiredParam<Real>("component", "component at which to compute this field ('0 for x, 1 for y, 2 for z')");
  return params;
}


MieFieldReals::MieFieldReals(const InputParameters & parameters) :
  AuxKernel(parameters),
  _a(getParam<Real>("a")),
  _omega(getParam<Real>("omega")),
  _c(getParam<Real>("c")),
  _epsilonI(getParam<Real>("epsilonI")),
  _sigmaI(getParam<Real>("sigmaI")),
  _epsilonII(getParam<Real>("epsilonII")),
  _sigmaII(getParam<Real>("sigmaII")),
  _L(getParam<Real>("L")),
  _order(getParam<Real>("order")),
  _component(getParam<Real>("component"))
{
}

Real
MieFieldReals::computeValue()
{
#ifdef FERRET_HAVE_BOOST_MATH_SPECIAL_FUNCTIONS
  Real x = _q_point[_qp](0);
  Real y = _q_point[_qp](1);
  Real z = _q_point[_qp](2);

  Complex k1I = Complex(0, 1) * (_omega / _c) * (_epsilonI + Complex(0, 1) * 4.0 * libMesh::pi * _sigmaI / _omega);
  Complex k2I = Complex(0, 1) * (_omega / _c);

  Complex ksqI = - k1I * k2I;

  Complex k1II = Complex(0, 1) * (_omega / _c) * (_epsilonII + Complex(0, 1) * 4.0 * libMesh::pi * _sigmaII / _omega);
  Complex k2II = Complex(0, 1) * (_omega / _c);

  Complex ksqII = - k1II * k2II;

  Real r =  std::pow(x, 2.0) + std::pow(y, 2.0) + std::pow(z, 2.0);
  Real th = std::acos(z / r);
  Real phi = std::atan(y / x); //this checks out
  Real q = 2 * libMesh::pi * _a / _L;

  Complex Ersj = 0.0;
  Complex Ethsj = 0.0;
  Complex Ephisj = 0.0;

  Complex Ex = 0.0;
  Complex Ey = 0.0;
  Complex Ez = 0.0;


  for(unsigned int j = 1; j < _order; ++j)
    {
      Complex eBj = ( (2.0 * j + 1.0) / (j *(j + 1.0 )) ) * Complex(0,1) * std::pow(Complex(0,1), j) * 1.0;
      //(k2I * std::sqrt(ksqII) * 1.0 - k2II * std::sqrt(ksqI) * 1.0) / (k2I * std::sqrt(ksqII) * 1.0 - k2II * std::sqrt(ksqI) * 1.0);

      //^^this line isn't working right... should switch to N and Q formulation... K will be useful if we want to solve for Mie effects on FE spheres.

      Complex mBj = ( (2.0 * j + 1.0) / (j *(j + 1.0 )) ) * Complex(0,1) * std::pow(Complex(0,1), j) * 1.0;
      //(k2I * std::sqrt(ksqII) * 1.0 - k2II * std::sqrt(ksqI) * 1.0)/(k2I * std::sqrt(ksqII) * 1.0 - k2II * std::sqrt(ksqI) * 1.0);

      Ersj += (std::pow(_L, 2.0) * std::cos(phi) / (4 * std::pow(libMesh::pi, 2.0) * r * r)) * j * (j + 1.0) * eBj * boost::math::legendre_p(j, 1, std::cos(th));

      Ethsj += (std::pow(_L, 2.0) * std::cos(phi) / (4 * std::pow(libMesh::pi, 2.0) * r)) * (
                                                                                             eBj * ((std::sin(th) * boost::math::legendre_p(j, 2, std::cos(th)) + std::cos(th) * boost::math::legendre_p(j, 1, std::cos(th)) )/ (std::sin(th) * std::sin(th))) * std::sin(th)
                                                                                             - Complex(0,1) * mBj * boost::math::legendre_p(j, 1, std::cos(th)) / std::sin(th) );

      Ephisj += (std::pow(_L, 2.0) * std::sin(phi) / (4 * std::pow(libMesh::pi, 2.0) * r)) * (
                                                                                              eBj / std::sin(th) * boost::math::legendre_p(j, 1, std::cos(th))
                                                                                              + Complex(0,1) * mBj * ((std::sin(th) * boost::math::legendre_p(j, 2, std::cos(th)) + std::cos(th) * boost::math::legendre_p(j, 1, std::cos(th)) )/ (std::sin(th) * std::sin(th))) * std::sin(th)
                                                                                              );
    }

  if (_component == 0)
    {
      Ex += Ethsj * std::cos(th) * std::cos(phi) +  Ersj * std::cos(phi) * std::sin(th) - Ephisj * std::sin(phi);
      return r * Ex.real();
    }
  else if (_component == 1)
    {
      Ey += Ersj * std::sin(th) * std::sin(phi) +  std::cos(phi) *(Ephisj + Ethsj * std::sin(phi) );
      return r * Ey.real();
    }
  else
    {
      Ez += Ersj * std::cos(th) - Ethsj * std::sin(phi);
      return r * Ez.real();
    }
#else
  return -1.0;
#endif
}
