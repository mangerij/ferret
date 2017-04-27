#include "ComputeElectroopticTensorBase.h"

template<>
InputParameters validParams<ComputeElectroopticTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeElectroopticTensorBase::ComputeElectroopticTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _electrooptic_tensor_name(_base_name + "electrooptic_tensor"),
   _electrooptic_tensor(declareProperty<RankThreeTensor>(_electrooptic_tensor_name))
{
}

void
ComputeElectroopticTensorBase::computeQpProperties()
{
  computeQpElectroopticTensor();
}
