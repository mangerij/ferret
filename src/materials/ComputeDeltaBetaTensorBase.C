
#include "ComputeDeltaBetaTensorBase.h"

template<>
InputParameters validParams<ComputeDeltaBetaTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeDeltaBetaTensorBase::ComputeDeltaBetaTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _delta_beta_tensor_name(_base_name + "delta_beta_tensor"),
   _delta_beta_tensor(declareProperty<RankTwoTensor>(_delta_beta_tensor_name)),
   _delta_beta_tensor_test_name(_base_name + "delta_beta_tensor_test"),
   _delta_beta_tensor_test(declareProperty<RankTwoTensor>(_delta_beta_tensor_test_name))
{
}

void
ComputeDeltaBetaTensorBase::computeQpProperties()
{
  computeQpDeltaBetaTensor();
}
