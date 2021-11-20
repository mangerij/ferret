#include "ComputeThermalConductivityTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeThermalConductivityTensor);
InputParameters
ComputeThermalConductivityTensor::validParams()
{
  InputParameters params = ComputeRotatedThermalConductivityTensorBase::validParams();
  params.addClassDescription("Compute a ThermalConductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("k_ij", "ThermalConductivity tensor for material");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeThermalConductivityTensor::ComputeThermalConductivityTensor(
    const InputParameters & parameters) : ComputeRotatedThermalConductivityTensorBase(parameters)
    // _kij(this->template getParam<std::vector<Real> >("k_ij"), (RankTwoTensor::FillMethod)(int)this->template getParam<MooseEnum>("fill_method"))
{
  // _kij(getParam<std::vector<Real> >("k_ij"), (RankTwoTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"));
  _kij.fillFromInputVector(getParam<std::vector<Real>>("k_ij"));
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate thermal conductivity tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _kij.rotate(R);
}
void
ComputeThermalConductivityTensor::computeQpThermalConductivityTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _thC_tensor[_qp] = _kij;
}
