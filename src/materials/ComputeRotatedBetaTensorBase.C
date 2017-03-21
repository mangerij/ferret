/*
 * @file   ComputeRotatedBetaTensor.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate photoelastic change to the refractive index
 *
 * where \Delta (1/n^2) = \Delta B_{ij} = p_{ijkl} \varepsilon_{kl}
 *
 * Note that B_{ij} = \epsilon_{ij}^{-1} in the principle axis frame.
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#include "ComputeRotatedBetaTensorBase.h"

template<>
InputParameters validParams<ComputeRotatedBetaTensorBase>()
{
  InputParameters params = validParams<ComputeBetaTensorBase>();
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}

ComputeRotatedBetaTensorBase::ComputeRotatedBetaTensorBase(const InputParameters & parameters) :
    ComputeBetaTensorBase(parameters),
    _Euler_angles(getParam<Real>("euler_angle_1"),
                  getParam<Real>("euler_angle_2"),
                  getParam<Real>("euler_angle_3"))
{
}

