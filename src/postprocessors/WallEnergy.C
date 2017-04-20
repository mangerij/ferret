/**
 * @file   WallEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <john.mangeri@uconn.edu
 * @date
 *
 * @brief This is a energy postprocessor that tracks the local
 *        gradients of the polar vector contributions to the energy.
 *
 */

#include "WallEnergy.h"

template<>
InputParameters validParams<WallEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
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

WallEnergy::WallEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11/G110")*_G110),
  _G12(getParam<Real>("G12/G110")*_G110),
  _G44(getParam<Real>("G44/G110")*_G110),
  _G44P(getParam<Real>("G44P/G110")*_G110),
  _len_scale(getParam<Real>("len_scale"))
{}

Real
WallEnergy::computeQpIntegral()
{
  return (0.5*_G11*(pow(_polar_x_grad[_qp](0),2)+pow(_polar_y_grad[_qp](1),2)+pow(_polar_z_grad[_qp](2),2))+
    _G12*(_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)+_polar_y_grad[_qp](1)*_polar_z_grad[_qp](2)+_polar_x_grad[_qp](0)*_polar_z_grad[_qp](2))+
    0.5*_G44*(pow(_polar_x_grad[_qp](1)+_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)+_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)+_polar_z_grad[_qp](0),2))+
	  0.5*_G44P*(pow(_polar_x_grad[_qp](1)-_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)-_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)-_polar_z_grad[_qp](0),2)))*_len_scale;
}
