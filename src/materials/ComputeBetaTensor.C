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
  params.addRequiredParam<Real>("n_o", "ordinary refractive index");
  params.addRequiredParam<Real>("n_e", "extraordinary refractive index");
  //params.addRequiredParam<std::vector<Real> >("b_ij", "impermeability tensor for material"); //diagonal3 for now
  //params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "diagonal3", "The fill method"); 
  return params;
}

ComputeBetaTensor::ComputeBetaTensor(const InputParameters & parameters) :
    ComputeRotatedBetaTensorBase(parameters),
   _no(getParam<Real>("n_o")),
   _ne(getParam<Real>("n_e"))
    //_bij(getParam<std::vector<Real> >("b_ij"), (RankTwoTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")) TODO: fix, how to fill_method a tensor?
{
}

void
ComputeBetaTensor::computeQpBetaTensor()
{
  //define the principle axis impermeability
  _beta_tensor[_qp](0,0) = 1.0 / (_no * _no * _no);
  _beta_tensor[_qp](0,1) = 0.0;
  _beta_tensor[_qp](0,2) = 0.0;
  _beta_tensor[_qp](1,0) = 0.0;
  _beta_tensor[_qp](1,1) = 1.0 / (_no * _no * _no);
  _beta_tensor[_qp](1,2) = 0.0;
  _beta_tensor[_qp](2,0) = 0.0;
  _beta_tensor[_qp](2,1) = 0.0;
  _beta_tensor[_qp](2,2) = 1.0 / (_ne * _ne * _ne);

  // Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles);
  _beta_tensor[_qp].rotate(R);
}


