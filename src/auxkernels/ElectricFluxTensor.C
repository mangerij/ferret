#include "ElectricFluxTensor.h"
#include "Material.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ElectricFluxTensor);

template<>
InputParameters validParams<ElectricFluxTensor> ()
{
InputParameters params = validParams<AuxKernel>();
params.addClassDescription("Electric flux generated");
params.addRequiredCoupledVar("T", "temperature");
params.addRequiredCoupledVar("potential_E_int", "electric potential");
params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
return params;
}

ElectricFluxTensor::ElectricFluxTensor(const InputParameters & parameters)
  : AuxKernel(parameters),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _ecC_tensor(getMaterialProperty<RankTwoTensor>("ecC_tensor")),
    _sbC_tensor(getMaterialProperty<RankTwoTensor>("sbC_tensor")),
    _component(getParam<unsigned int>("component"))
{

}

Real
ElectricFluxTensor::computeValue()
//tensor inclusion
{
  Real sum = 0.0;
  for (unsigned int i = 0, j = 0, k = 0; i < 3 && j < 3 && k < 3; ++i, ++j, ++k)
  {
  // sum += -_ecC_tensor[_qp](i,0) * _potential_E_int_grad[_qp](0) -
  //         _sbC_tensor[_qp](i,k) * _ecC_tensor[_qp](k,0) * _T_grad[_qp](0) -
  //         _ecC_tensor[_qp](i,1) * _potential_E_int_grad[_qp](1) -
  //         _sbC_tensor[_qp](i,k) * _ecC_tensor[_qp](k,1) * _T_grad[_qp](1) -
  //         _ecC_tensor[_qp](i,2) * _potential_E_int_grad[_qp](2) -
  //         _sbC_tensor[_qp](i,k) * _ecC_tensor[_qp](k,2) * _T_grad[_qp](2);

     sum += -_ecC_tensor[_qp](i,j) * _potential_E_int_grad[_qp](_component) -
             _sbC_tensor[_qp](i,j) * _ecC_tensor[_qp](j,k) * _T_grad[_qp](_component);
}
return sum;
}
