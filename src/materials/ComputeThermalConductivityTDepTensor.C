#include "ComputeThermalConductivityTDepTensor.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeThermalConductivityTDepTensor);
InputParameters
ComputeThermalConductivityTDepTensor::validParams()
{
  InputParameters params = ComputeRotatedThermalConductivityTensorBase::validParams();
  params.addClassDescription("Compute a ThermalConductivity tensor.");
  params.addRequiredParam<std::vector<Real>>("ak_ij", "ThermalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("bk_ij", "ThermalConductivity tensor for material");
  params.addRequiredParam<std::vector<Real>>("ck_ij", "ThermalConductivity tensor for material");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<MooseEnum>("fill_method", RankTwoTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}
ComputeThermalConductivityTDepTensor::ComputeThermalConductivityTDepTensor(
    const InputParameters & parameters) : ComputeRotatedThermalConductivityTensorBase(parameters),
   _T(coupledValue("T"))
{
  _akij.fillFromInputVector(getParam<std::vector<Real>>("ak_ij"));
  _bkij.fillFromInputVector(getParam<std::vector<Real>>("bk_ij"));
  _ckij.fillFromInputVector(getParam<std::vector<Real>>("ck_ij"));

  RotationTensor R(_Euler_angles);
  _akij.rotate(R);
  _bkij.rotate(R);
  _ckij.rotate(R);
}
void
ComputeThermalConductivityTDepTensor::computeQpThermalConductivityTensor()
{
  _thC_tensor[_qp] = _akij + _bkij*_T[_qp] + _ckij*_T[_qp]*_T[_qp];
}
