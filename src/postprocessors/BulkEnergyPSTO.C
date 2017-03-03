/**
 * @file   BulkEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun  4 15:01:35 2013
 *
 * @brief
 *
 *
 */


#include "BulkEnergyPSTO.h"

template<>
InputParameters validParams<BulkEnergyPSTO>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion"); 
  params.addRequiredParam<Real>("alpha2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha3", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha4", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha5", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha6", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

BulkEnergyPSTO::BulkEnergyPSTO(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha2(getParam<Real>("alpha2")),
  _alpha3(getParam<Real>("alpha3")),
  _alpha4(getParam<Real>("alpha4")),
  _alpha5(getParam<Real>("alpha5")),
  _x1(getParam<Real>("x1")),
  _x2(getParam<Real>("x2")),
  _x3(getParam<Real>("x3")),
  _x4(getParam<Real>("x4")),
  _x5(getParam<Real>("x5")),
  _x6(getParam<Real>("x6")),
  _epsilon(getParam<Real>("epsilon")),
  _T(getParam<Real>("T")),
  _Tc(getParam<Real>("Tc"))

{
}

Real
BulkEnergyPSTO::computeQpIntegral()
{
  return 
  _alpha1 * (_T-_Tc) * (std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0)) + _alpha2 * (std::pow(_polar_x[_qp], 4.0) +std::pow(_polar_y[_qp], 4.0)) + _alpha3 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + _alpha4 * (std::pow(_polar_x[_qp], 6.0) + std::pow(_polar_y[_qp], 6.0)) + _alpha5 * (std::pow(_polar_x[_qp], 4.0) * std::pow(_polar_y[_qp], 2.0) + std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 4.0)) + (_x1 * (std::pow(_polar_x[_qp], 2.0) +std::pow(_polar_y[_qp], 2.0)) + _x2 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) + _x3 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 4.0)) * _epsilon + (_x4 * (std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0)) + _x5 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0))+ _x6 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0)) * _epsilon * _epsilon; 

}
