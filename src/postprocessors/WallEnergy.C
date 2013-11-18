/**
 * @file   WallEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jul 30 16:35:05 2013
 *
 * @brief
 *
 *
 */

#include "WallEnergy.h"

template<>
InputParameters validParams<WallEnergy>()
{
  //TODO: inherit from an appropriate postprocessor
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<Real>("G110"," "); //FIXME: Give me an explanation
  params.addRequiredParam<Real>("G11/G110"," ");
  params.addRequiredParam<Real>("G12/G110"," ");
  params.addRequiredParam<Real>("G44/G110"," ");
  params.addRequiredParam<Real>("G44P/G110"," ");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  params.addParam<Real>("energy_scale",1.0,"energy scale");
  return params;
}

WallEnergy::WallEnergy(const std::string & name, InputParameters parameters) :
  ElementIntegralPostprocessor(name, parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11/G110")*_G110),
  _G12(getParam<Real>("G12/G110")*_G110),
  _G44(getParam<Real>("G44/G110")*_G110),
  _G44P(getParam<Real>("G44P/G110")*_G110),
  _len_scale(getParam<Real>("len_scale")),
  _energy_scale(getParam<Real>("energy_scale"))
{}

Real
WallEnergy::computeQpIntegral()
{
  return (0.5*_G11*(pow(_polar_x_grad[_qp](0),2)+pow(_polar_y_grad[_qp](1),2)+pow(_polar_z_grad[_qp](2),2))+
    _G12*(_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)+_polar_y_grad[_qp](1)*_polar_z_grad[_qp](2)+_polar_x_grad[_qp](0)*_polar_z_grad[_qp](2))+
    0.5*_G44*(pow(_polar_x_grad[_qp](1)+_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)+_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)+_polar_z_grad[_qp](0),2))+
	  0.5*_G44P*(pow(_polar_x_grad[_qp](1)-_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)-_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)-_polar_z_grad[_qp](0),2)))*_len_scale*_energy_scale;
}
