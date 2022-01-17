#include "ComputeElectricalConductivityTDepTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeElectricalConductivityTDepTensor);
InputParameters
ComputeElectricalConductivityTDepTensor::validParams()
{
  InputParameters params = ComputeRotatedElectricalConductivityTensorBase::validParams();
  params.addClassDescription("Compute a ElectricalConductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("asg_ij", "ElectricalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("bsg_ij", "ElectricalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("csg_ij", "ElectricalConductivity tensor for material");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeElectricalConductivityTDepTensor::ComputeElectricalConductivityTDepTensor(
    const InputParameters & parameters)
  : ComputeRotatedElectricalConductivityTensorBase(parameters),
   _T(coupledValue("T"))
{
  _asgij.fillFromInputVector(getParam<std::vector<Real>>("asg_ij"));
  _bsgij.fillFromInputVector(getParam<std::vector<Real>>("bsg_ij"));
  _csgij.fillFromInputVector(getParam<std::vector<Real>>("csg_ij"));

  RotationTensor R(_Euler_angles);

  _asgij.rotate(R);
  _bsgij.rotate(R);
  _csgij.rotate(R);
}

void
ComputeElectricalConductivityTDepTensor::computeQpElectricalConductivityTensor()
{
  ///Assign an electrical conductivity tensor at a given quad point. This will be reworked eventually for constant _qp.
  _ecC_tensor[_qp] = _asgij+_bsgij*_T[_qp] + _csgij*_T[_qp]*_T[_qp];
}
