/**
 * @file   ComputeBetaTensor.C
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

#include "ComputeBetaTensor.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

template<>
InputParameters validParams<ComputeBetaTensor>()
{
  InputParameters params = validParams<ComputeRotatedBetaTensorBase>();
  params.addClassDescription("Compute the impermeability tensor.");
  params.addRequiredParam<Real>("n_a", "alpha refractive index");
  params.addRequiredParam<Real>("n_b", "beta refractive index");
  params.addRequiredParam<Real>("n_g", "gamma refractive index");
  //params.addRequiredParam<std::vector<Real> >("b_ij", "impermeability tensor for material"); //diagonal3 for now
  //params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "diagonal3", "The fill method"); 
  return params;
}

ComputeBetaTensor::ComputeBetaTensor(const InputParameters & parameters) :
    ComputeRotatedBetaTensorBase(parameters),
   _na(getParam<Real>("n_a")),
   _nb(getParam<Real>("n_b")),
   _ng(getParam<Real>("n_g"))
    //_bij(getParam<std::vector<Real> >("b_ij"), (RankTwoTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")) TODO: fix, how to fill_method a tensor?
{
}

void
ComputeBetaTensor::computeQpBetaTensor()
{
  // Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles);
  //define the principle axis impermeability rotated
  _beta_tensor[_qp](0,0) = R(0,0) * R(0,0) / (_na * _na) + R(1,0) * R(1,0) / (_nb * _nb) + R(2,0) * R(2,0) / (_ng * _ng);//
  _beta_tensor[_qp](0,1) = R(0,0) * R(0,1) / (_na * _na) + R(1,0) * R(1,1) / (_nb * _nb) + R(2,1) * R(2,0) / (_ng * _ng);//
  _beta_tensor[_qp](0,2) = R(0,0) * R(0,2) / (_na * _na) + R(1,0) * R(1,2) / (_nb * _nb) + R(2,0) * R(2,2) / (_ng * _ng);//
  _beta_tensor[_qp](1,0) = R(0,1) * R(0,0) / (_na * _na) + R(1,1) * R(1,0) / (_nb * _nb) + R(2,0) * R(2,1) / (_ng * _ng);//
  _beta_tensor[_qp](1,1) = R(0,1) * R(0,1) / (_na * _na) + R(1,1) * R(1,1) / (_nb * _nb) + R(2,1) * R(2,1) / (_ng * _ng);//
  _beta_tensor[_qp](1,2) = R(0,1) * R(0,2) / (_na * _na) + R(1,1) * R(1,2) / (_nb * _nb) + R(2,1) * R(2,2) / (_ng * _ng);//
  _beta_tensor[_qp](2,0) = R(0,2) * R(0,0) / (_na * _na) + R(1,0) * R(1,2) / (_nb * _nb) + R(2,0) * R(2,2) / (_ng * _ng);//
  _beta_tensor[_qp](2,1) = R(0,1) * R(0,2) / (_na * _na) + R(1,1) * R(1,2) / (_nb * _nb) + R(2,1) * R(2,2) / (_ng * _ng);//
  _beta_tensor[_qp](2,2) = R(0,2) * R(0,2) / (_na * _na) + R(1,2) * R(1,2) / (_nb * _nb) + R(2,2) * R(2,2) / (_ng * _ng);//

  //Moose::out << "\n B"; std::cout << 0; std::cout << 0; Moose::out << " = "; std::cout << _beta_tensor[_qp](0,0);
}


