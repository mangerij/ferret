
#include "ComputeBetaTensorBase.h"

template<>
InputParameters validParams<ComputeBetaTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeBetaTensorBase::ComputeBetaTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _beta_tensor_name(_base_name + "beta_tensor"),
   _beta_tensor(declareProperty<RankTwoTensor>(_beta_tensor_name))
{
}

void
ComputeBetaTensorBase::computeQpProperties()
{
  computeQpBetaTensor();
}
