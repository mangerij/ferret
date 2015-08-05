/**
 * @file   BulkEnergyFourth.C
 * @author J. Mangeri <mangerij@anl.gov>
 * @date   Tue Aug  215:01:35 2015
 *
 * @brief
 *
 *
 */


#include "BulkEnergyFourth.h"

template<>
InputParameters validParams<BulkEnergyFourth>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1","The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha11","The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha12","The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

BulkEnergyFourth::BulkEnergyFourth(const std::string & name, InputParameters parameters) :
  ElementIntegralPostprocessor(name, parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha12(getParam<Real>("alpha12")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkEnergyFourth::computeQpIntegral()
{
  return ((pow(_polar_x[_qp],2)+pow(_polar_y[_qp],2)+pow(_polar_z[_qp],2))+
    _alpha12*(1/_alpha1)*(pow(_polar_x[_qp],2)*pow(_polar_y[_qp],2)+
	      pow(_polar_y[_qp],2)*pow(_polar_z[_qp],2)+
	      pow(_polar_x[_qp],2)*pow(_polar_z[_qp],2))
	  )*pow(_len_scale,3)*_alpha1;
}
