#include "ComputeElasticityTensor.h"
#include "ComputeElectrostrictiveTensorBase.h"

template<>
InputParameters validParams<ComputeElectrostrictiveTensorBase>()
{
  InputParameters params = validParams<ComputeElasticityTensor>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeElectrostrictiveTensorBase::ComputeElectrostrictiveTensorBase(const InputParameters & parameters) :
    ComputeElasticityTensor(parameters)
//    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
//    _electrostrictive_tensor_name(_base_name + "electrostrictive_tensor"),
//    _electrostrictive_tensor(declareProperty<ElectrostrictiveTensorR4>(_electrostrictive_tensor_name))

{
}

void
ComputeElectrostrictiveTensorBase::computeQpProperties()
{
  computeQpElectrostrictiveTensor();
}
