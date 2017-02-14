#include "ComputeElasticityTensor.h"
#include "ComputePhotostrictiveTensorBase.h"

template<>
InputParameters validParams<ComputePhotostrictiveTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputePhotostrictiveTensorBase::ComputePhotostrictiveTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _photostrictive_tensor_name(_base_name + "photostrictive_tensor"),
   _photostrictive_tensorP_name(_base_name + "photostrictive_tensorQ"),
   _photostrictive_tensor(declareProperty<RankFourTensor>(_photostrictive_tensor_name)),
   _photostrictive_tensorP(declareProperty<RankFourTensor>(_photostrictive_tensorP_name))
{
}

void
ComputePhotostrictiveTensorBase::computeQpProperties()
{
  computeQpPhotostrictiveTensor();
}
