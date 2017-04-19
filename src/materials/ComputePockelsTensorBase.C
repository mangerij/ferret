#include "ComputePockelsTensorBase.h"

template<>
InputParameters validParams<ComputePockelsTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputePockelsTensorBase::ComputePockelsTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _pockels_tensor_name(_base_name + "pockels_tensor"),
   _pockels_tensor(declareProperty<RankFourTensor>(_pockels_tensor_name))
{
}

void
ComputePockelsTensorBase::computeQpProperties()
{
  computeQpPockelsTensor();
}
