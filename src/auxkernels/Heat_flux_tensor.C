#include "Heat_flux_tensor.h"
#include "Material.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", Heat_flux_tensor);

template <>
InputParameters
validParams<Heat_flux_tensor>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("heat flux generated");
  params.addRequiredCoupledVar("T", "temperature");
  params.addRequiredCoupledVar("potential_E_int", "electric potential");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction the variable "
                                        "this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

Heat_flux_tensor::Heat_flux_tensor(const InputParameters & parameters)
  : AuxKernel(parameters),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _thC_tensor(getMaterialProperty<RankTwoTensor>("thC_tensor")), // added for tensor application
    _ecC_tensor(getMaterialProperty<RankTwoTensor>("ecC_tensor")),
    _sbC_tensor(getMaterialProperty<RankTwoTensor>("sbC_tensor")),
    _component(getParam<unsigned int>("component")),
    _len_scale(getParam<Real>("len_scale"))
{
}

Real
Heat_flux_tensor::computeValue()
// tensor inclusion
{
  Real sum = 0.0;
  for (unsigned int i = 0, k = 0, l = 0; i < 3 && k < 3 && l < 3; ++i, ++k, ++l)
  {
    sum +=
        (-_thC_tensor[_qp](i, 0) * _T_grad[_qp](0) -
         _sbC_tensor[_qp](i, k) * _ecC_tensor[_qp](k, 0) * _T[_qp] * _potential_E_int_grad[_qp](0) -
         _sbC_tensor[_qp](i, k) * _sbC_tensor[_qp](k, l) * _ecC_tensor[_qp](l, 0) * _T[_qp] *
             _T_grad[_qp](0) -
         _thC_tensor[_qp](i, 1) * _T_grad[_qp](1) -
         _sbC_tensor[_qp](i, k) * _ecC_tensor[_qp](k, 1) * _T[_qp] * _potential_E_int_grad[_qp](1) -
         _sbC_tensor[_qp](i, k) * _sbC_tensor[_qp](k, l) * _ecC_tensor[_qp](l, 1) * _T[_qp] *
             _T_grad[_qp](1) -
         _thC_tensor[_qp](i, 2) * _T_grad[_qp](2) -
         _sbC_tensor[_qp](i, k) * _ecC_tensor[_qp](k, 2) * _T[_qp] * _potential_E_int_grad[_qp](2) -
         _sbC_tensor[_qp](i, k) * _sbC_tensor[_qp](k, l) * _ecC_tensor[_qp](l, 2) * _T[_qp] *
             _T_grad[_qp](2)) *
        _len_scale;
  }
  return sum;
}
