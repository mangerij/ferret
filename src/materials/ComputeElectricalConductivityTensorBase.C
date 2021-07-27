
#include "ComputeElectricalConductivityTensorBase.h"

template <>
InputParameters
validParams<ComputeElectricalConductivityTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple mechanics material systems on "
      "the same block, i.e. for multiple phases");
  return params;
}

ComputeElectricalConductivityTensorBase::ComputeElectricalConductivityTensorBase(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _ecC_tensor_name(_base_name + "ecC_tensor"),
    _ecC_tensor(declareProperty<RankTwoTensor>(_ecC_tensor_name))
{
}

void
ComputeElectricalConductivityTensorBase::computeQpProperties()
{
  computeQpElectricalConductivityTensor();
}
