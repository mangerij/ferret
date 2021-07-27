
#include "ComputeSeebeckTensorBase.h"

template <>
InputParameters
validParams<ComputeSeebeckTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple mechanics material systems on "
      "the same block, i.e. for multiple phases");
  return params;
}

ComputeSeebeckTensorBase::ComputeSeebeckTensorBase(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _sbC_tensor_name(_base_name + "sbC_tensor"),
    _sbC_tensor(declareProperty<RankTwoTensor>(_sbC_tensor_name))
{
}

void
ComputeSeebeckTensorBase::computeQpProperties()
{
  computeQpSeebeckTensor();
}
