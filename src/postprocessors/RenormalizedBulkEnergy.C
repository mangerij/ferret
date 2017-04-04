/**
 * @file   RenormalizedBulkEnergy.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief
 *
 *
 */


#include "RenormalizedBulkEnergy.h"

template<>
InputParameters validParams<RenormalizedBulkEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("epsilon",0.0,"the strain coupling");
  params.addParam<Real>("T",0.0,"the temperature");
  return params;
}

RenormalizedBulkEnergy::RenormalizedBulkEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _epsilon(getParam<Real>("epsilon")),
  _T(getParam<Real>("T"))
{
}

Real
RenormalizedBulkEnergy::computeQpIntegral()
{
  return 0.488524 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + 0.77016 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) - 3.7 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 0.132318 * std::pow(_polar_z[_qp], 4.0) + 0.253636 * (std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0)) + 0.26 * ( std::pow(_polar_x[_qp], 6.0) + std::pow(_polar_y[_qp], 6.0) + std::pow(_polar_z[_qp], 6.0)) + 0.61 * ((std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) *std::pow(_polar_z[_qp], 2.0) + (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0)) *std::pow(_polar_y[_qp], 2.0) + (std::pow(_polar_y[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0)) *std::pow(_polar_x[_qp], 2.0)) + (std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0)) * (-4.92223 * (-765.1 + _T) - 19.0909 * _epsilon) + std::pow(_polar_z[_qp], 2.0) * (-4.92223 * (-765.1 + _T) + 15.7576 * _epsilon);
  //0.488524 Px^2 Py^2 + 0.77016 (Px^4 + Py^4) - 3.7 Px^2 Py^2 Pz^2 +  0.132318 Pz^4 + 0.253636 (Px^2 Pz^2 + Py^2 Pz^2) +  0.26 (Px^6 + Py^6 + Pz^6) +  0.61 ((Px^4 + Py^4) Pz^2 + Py^2 (Px^4 + Pz^4) +     Px^2 (Py^4 + Pz^4)) + (Px^2 + Py^2) (0.000229917 -     19.0909 \[Epsilon]) + 303.03 \[Epsilon]^2 +  Pz^2 (0.000229917 + 15.7576 \[Epsilon])
}
