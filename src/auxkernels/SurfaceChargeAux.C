#include "SurfaceChargeAux.h"

template<>

InputParameters validParams<SurfaceChargeAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}


SurfaceChargeAux::SurfaceChargeAux( const std::string & name, InputParameters parameters ) :
  AuxKernel( name, parameters ),
  _normals(_var.normals()),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{
}

Real
SurfaceChargeAux::computeValue()
{
  RealVectorValue P(_polar_x[_qp],_polar_y[_qp],_polar_z[_qp]);
  RealVectorValue n(_normals[_qp]);
    return P*n;
}


