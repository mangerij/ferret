/**
 * @file   BulkEnergyCoupledT.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @brief
 *
 *
 */


#include "BulkEnergyCoupledT.h"

template<>
InputParameters validParams<BulkEnergyCoupledT>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("temperature", "The temperature at the grid point");
  params.addRequiredParam<Real>("alpha0", "The coefficients of the Landau expansion"); 
  params.addRequiredParam<Real>("alpha11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("Tc", "Transition temperature of unstrained ferroelectric");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

BulkEnergyCoupledT::BulkEnergyCoupledT(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _temperature(coupledValue("temperature")),
  _alpha0(getParam<Real>("alpha0")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha12(getParam<Real>("alpha12")),
  _alpha111(getParam<Real>("alpha111")),
  _alpha112(getParam<Real>("alpha112")),
  _alpha123(getParam<Real>("alpha123")),
  _Tc(getParam<Real>("Tc")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkEnergyCoupledT::computeQpIntegral()
{
  return (
    _alpha0 * (_temperature[_qp] - _Tc) * (std::pow(_polar_x[_qp], 2) + std::pow(_polar_y[_qp], 2) + std::pow(_polar_z[_qp] ,2))
  + _alpha11 * (std::pow(_polar_x[_qp], 4) + std::pow(_polar_y[_qp], 4) + std::pow(_polar_z[_qp], 4))
    + _alpha12 * (std::pow(_polar_x[_qp], 2) * std::pow(_polar_y[_qp], 2)+
	      std::pow(_polar_y[_qp], 2) * std::pow(_polar_z[_qp], 2)+
	      std::pow(_polar_x[_qp], 2) * std::pow(_polar_z[_qp], 2))+
    _alpha111 * (std::pow(_polar_x[_qp], 6) + std::pow(_polar_y[_qp], 6) + std::pow(_polar_z[_qp], 6))+
    _alpha112 * (std::pow(_polar_x[_qp], 4) * (std::pow(_polar_y[_qp], 2) + std::pow(_polar_z[_qp], 2))
	      + std::pow(_polar_y[_qp], 4) * (std::pow(_polar_z[_qp], 2) + std::pow(_polar_x[_qp], 2))
	      + std::pow(_polar_z[_qp], 4) * (std::pow(_polar_x[_qp], 2) + std::pow(_polar_y[_qp], 2)))+
	  _alpha123 * (pow(_polar_x[_qp], 2) * std::pow(_polar_y[_qp], 2) * std::pow(_polar_z[_qp], 2))) * std::pow(_len_scale,3);
}
