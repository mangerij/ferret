/****************************************************************/
/* Computes a rank 4 electrostrictive tensor                    */
/****************************************************************/

#include "ComputeElasticityTensor.h"
#include "ComputeRotatedElasticityTensorBase.h"
#include "ComputePhotostrictiveTensor.h"
#include "RotationTensor.h"
#include "RankFourTensor.h"

template<>
InputParameters validParams<ComputePhotostrictiveTensor>()
{
  InputParameters params = validParams<ComputeRotatedPhotostrictiveTensorBase>();
  params.addClassDescription("Compute a photostrictive tensor.");
  params.addRequiredParam<std::vector<Real> >("P_mnkl", "elasto-optic tensor for material");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputePhotostrictiveTensor::ComputePhotostrictiveTensor(const InputParameters & parameters) :
    ComputeRotatedPhotostrictiveTensorBase(parameters),
    _Pmnkl(getParam<std::vector<Real> >("P_mnkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _Pmnkl.rotate(R);
}

void
ComputePhotostrictiveTensor::computeQpPhotostrictiveTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _photostrictive_tensor[_qp] = _Pmnkl;
}


//void
//ComputePhotostrictiveTensor::computeQpUnstrainedRefractiveIndex()
//{
//  // Assume that n_e is along the z-axis for now
//  // note the regular birefringence is quantified by _ne - _no
//  RealVectorValue n(_no, _no, _ne); 

//  // Rotate the indicatrix such that it is aligned with the crystallographic direction A_i = R_{ij} A_j
//  RealVectorValue nR(R(0, 0) * n(0) + R(0, 1) * n(1) + R(0, 2) * n(2), R(1, 0) * n(0) + R(1, 1) * n(1) + R(1, 2) * n(2), R(2, 0) * n(0) + R(2, 1) * n(1) + R(2, 2) * n(2));1
//}

//void
//ComputePhotostrictiveTensor::computeQpStrainedRefractiveIndex()
//{
//  // Assume that n_e is along the z-axis for now
//  // note the regular birefringence is quantified by _ne - _no
//  RealVectorValue n(_no, _no, _ne); 
//
//  // Rotate the indicatrix such that it is aligned with the crystallographic direction A_i = R_{ij} A_j
//  RealVectorValue nR(R(0, 0) * n(0) + R(0, 1) * n(1) + R(0, 2) * n(2), R(1, 0) * n(0) + R(1, 1) * n(1) + R(1, 2) * n(2), R(2, 0) * n(0) + R(2, 1) * n(1) + R(2, 2) * n(2));1
//
//
//  then store in an aux kernel 
//}
