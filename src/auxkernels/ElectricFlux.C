#include "ElectricFlux.h"
#include "Material.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ElectricFlux);

template <>
InputParameters
validParams<ElectricFlux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Electric potential generated due to heat flux");
  params.addRequiredCoupledVar("T", "temperature");
  params.addRequiredCoupledVar("potential_E_int", "electric potential");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction the variable "
                                        "this kernel acts in. (0 for x, 1 for y, 2 for z)");
  return params;
}

ElectricFlux::ElectricFlux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _ecC(getMaterialProperty<Real>("ecC")),
    _sbC(getMaterialProperty<Real>("sbC")),
    _component(getParam<unsigned int>("component"))
{
}

Real
ElectricFlux::computeValue()
{
  return -_ecC[_qp] * _potential_E_int_grad[_qp](_component) -
         _sbC[_qp] * _ecC[_qp] * _T_grad[_qp](_component);
}
