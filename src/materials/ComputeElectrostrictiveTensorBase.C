/**
 * @file   ComputeElectrostrictiveTensorBase.C
 * @author J. Mangeri <john.mangeri@uconn.edu
 *
 * @brief  Base class for electrostrictive material
 *
 */


#include "ComputeElasticityTensor.h"
#include "ComputeElectrostrictiveTensorBase.h"

template<>
InputParameters validParams<ComputeElectrostrictiveTensorBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeElectrostrictiveTensorBase::ComputeElectrostrictiveTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _electrostrictive_tensor_name(_base_name + "electrostrictive_tensor"),
   _electrostrictive_coefficients_name(_base_name + "electrostrictive_coefficients"),
   _electrostrictive_tensor(declareProperty<RankFourTensor>(_electrostrictive_tensor_name)),
   _electrostrictive_coefficients(declareProperty<RankFourTensor>(_electrostrictive_coefficients_name))
{
}

void
ComputeElectrostrictiveTensorBase::computeQpProperties()
{
  computeQpElectrostrictiveTensor();
}
