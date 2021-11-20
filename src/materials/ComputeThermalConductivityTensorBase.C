
#include "ComputeThermalConductivityTensorBase.h"

InputParameters
ComputeThermalConductivityTensorBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple mechanics material systems on "
      "the same block, i.e. for multiple phases");
  return params;
}

ComputeThermalConductivityTensorBase::ComputeThermalConductivityTensorBase(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _thC_tensor_name(_base_name + "thC_tensor"),
    _thC_tensor(declareProperty<RankTwoTensor>(_thC_tensor_name))
{
}

void
ComputeThermalConductivityTensorBase::computeQpProperties()
{
  computeQpThermalConductivityTensor();
}
