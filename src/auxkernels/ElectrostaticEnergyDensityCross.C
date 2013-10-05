/**
 * @file   ElectrostaticEnergyDensityCross.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 17:14:17 2013
 *
 * @brief
 *
 *
 */


#include "ElectrostaticEnergyDensityCross.h"

template<>
InputParameters validParams<ElectrostaticEnergyDensityCross>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential", "The electrostatic potential");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

ElectrostaticEnergyDensityCross::ElectrostaticEnergyDensityCross(const std::string & name, InputParameters parameters) :
  AuxKernel(name, parameters),
  _potential_grad(coupledGradient("potential")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{}


Real
ElectrostaticEnergyDensityCross::computeValue()
{
  // Real r;
  // Real a[3],b[3];
  // r=_potential_grad[_qp](0)*_polar_x[_qp]+ _potential_grad[_qp](1)*_polar_y[_qp]+ _potential_grad[_qp](2)*_polar_z[_qp];
  // std::cout<<"r="<<r<<"    ("<<_potential_grad[_qp](0)<<","<<_potential_grad[_qp](1)<<","<<_potential_grad[_qp](2)<<")   "<<"("<<_polar_x[_qp]<<","<<_polar_y[_qp]<<","<<_polar_z[_qp]<<")\n";
  // a[0]=_potential_grad[_qp](0);a[1]=_potential_grad[_qp](1);a[2]=_potential_grad[_qp](2);
  // b[0]=_polar_x[_qp];b[1]=_polar_y[_qp];b[2]=_polar_z[_qp];

  // r=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  // std::cout<<"r="<<r<<"    ("<<a[0]<<","<<a[1]<<","<<a[2]<<")   "<<"("<<b[0]<<","<<b[1]<<","<<b[2]<<")\n\n";
  return _potential_grad[_qp](0)*_polar_x[_qp]+ _potential_grad[_qp](1)*_polar_y[_qp]+ _potential_grad[_qp](2)*_polar_z[_qp];
}
