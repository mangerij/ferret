
#include "ComputePolarOpticTensorBase.h"

template<>
InputParameters validParams<ComputePolarOpticTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputePolarOpticTensorBase::ComputePolarOpticTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _delta_PO_tensor_name(_base_name + "delta_PO_tensor"),
   _delta_PO_tensor(declareProperty<RankTwoTensor>("delta_PO_tensor"))
{
}

void
ComputePolarOpticTensorBase::computeQpProperties()
{
  computeQpPolarOpticTensor();
}
