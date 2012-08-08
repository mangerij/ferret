#include "VortexSurfaceEnergy.h"

template<>
InputParameters validParams<VortexSurfaceEnergy>()
{
  InputParameters params = validParams<SideIntegral>();
  params.addParam<Real>("a", 1.0, "The location of the vortex along the x-axis.");
  params.addParam<Real>("verticality", 1e-8, "'Verticality' tolerance: if the z-component of the normal is > than verticality, it is not horizontal.");
  return params;
}

VortexSurfaceEnergy::VortexSurfaceEnergy(const std::string & name, InputParameters parameters) :
    SideIntegral(name, parameters),
    _a(parameters.get<Real>("a")),
    _verticality(parameters.get<Real>("verticality"))
{}

Real
VortexSurfaceEnergy::computeQpIntegral()
{
  Real  res     = 0.0;
  Point xypoint = _q_point[_qp];
  xypoint(2)    = 0.0;
  Real  r       = sqrt(xypoint*xypoint);
  Real  cphi    = xypoint(0)/r;
  Real  sphi    = xypoint(1)/r;
  // \sigma = 2a \sin{\phi} \frac{a \cos{\phi} - 1}{1 + a^2 - 2a\cos{\phi}}
  Real  sigma   = 2.0*_a*sphi*(_a*cphi-1.0)/(1.0+_a*_a-2.0*_a*cphi);
  if(fabs(_normals[_qp](2)) < _verticality)
    res = _u[_qp]*sigma;  
  return res;
}
