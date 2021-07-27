#include "ComputeElectricalConductivityTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeElectricalConductivityTensor);
template <>
InputParameters
validParams<ComputeElectricalConductivityTensor>()
{
  InputParameters params = validParams<ComputeRotatedElectricalConductivityTensorBase>();
  params.addClassDescription("Compute a ElectricalConductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("g_ij", "ElectricalConductivity tensor for material");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeElectricalConductivityTensor::ComputeElectricalConductivityTensor(
    const InputParameters & parameters)
  : ComputeRotatedElectricalConductivityTensorBase(parameters)
// _gamma_ij(getParam<std::vector<Real> >("g_ij"),
// (RankTwoTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  _gij.fillFromInputVector(getParam<std::vector<Real>>("g_ij"));
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _gij.rotate(R);
}
void
ComputeElectricalConductivityTensor::computeQpElectricalConductivityTensor()
{
  ///Assign an electrical conductivity tensor at a given quad point. This will be reworked eventually for constant _qp.
  _ecC_tensor[_qp] = _gij;
}
