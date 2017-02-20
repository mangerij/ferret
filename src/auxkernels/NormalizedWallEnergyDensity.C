/**
 * @file   NormalizedWallEnergyDensity.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief calculate the gradient energy density:
 *
 *
 */
#include "NormalizedWallEnergyDensity.h"
#include <math.h>

template<>
InputParameters validParams<NormalizedWallEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("G110","Domain wall penalty coefficients");
  params.addRequiredParam<Real>("G11/G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G12/G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G44/G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G44P/G110","Ratio of domain wall penalty coefficients");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

NormalizedWallEnergyDensity::NormalizedWallEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _polar_x_var(coupled("polar_x")),
  _polar_y_var(coupled("polar_y")),
  _polar_z_var(coupled("polar_z")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11/G110")*_G110),
  _G12(getParam<Real>("G12/G110")*_G110),
  _G44(getParam<Real>("G44/G110")*_G110),
  _G44P(getParam<Real>("G44P/G110")*_G110),
  _len_scale(getParam<Real>("len_scale"))
{}

Real
NormalizedWallEnergyDensity::computeValue()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  return (1/sqrt(w*w))*(0.5*_G11*(pow(_polar_x_grad[_qp](0),2)+pow(_polar_y_grad[_qp](1),2)+pow(_polar_z_grad[_qp](2),2))+
    _G12*(_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)+_polar_y_grad[_qp](1)*_polar_z_grad[_qp](2)+_polar_x_grad[_qp](0)*_polar_z_grad[_qp](2))+
    0.5*_G44*(pow(_polar_x_grad[_qp](1)+_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)+_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)+_polar_z_grad[_qp](0),2))+
    0.5*_G44P*(pow(_polar_x_grad[_qp](1)-_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)-_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)-_polar_z_grad[_qp](0),2)))*_len_scale;
}
