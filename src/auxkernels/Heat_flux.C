#include "Heat_flux.h"
#include "Material.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", Heat_flux);

template <>
InputParameters
validParams<Heat_flux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Electric potential generated due to heat flux");
  params.addRequiredCoupledVar("T", "temperature");
  params.addRequiredCoupledVar("potential_E_int", "electric potential");
  params.addParam<MaterialPropertyName>(
      "thC", "Thermal Conductivity", "Property name of the thermal conductivity material property");
  params.addParam<MaterialPropertyName>(
      "ecC",
      "Electrical Conductivity",
      "Property name of the electrical conductivity material property");
  params.addParam<MaterialPropertyName>(
      "sbC", "Seebeck coefficient", "Property name of the Seebeck coefficient material property");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction the variable "
                                        "this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

Heat_flux::Heat_flux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _thC(getMaterialProperty<Real>("thC")),
    _ecC(getMaterialProperty<Real>("ecC")),
    _sbC(getMaterialProperty<Real>("sbC")),
    _component(getParam<unsigned int>("component")),
    _len_scale(getParam<Real>("len_scale"))
{
}

Real
Heat_flux::computeValue()
{
  return (-_thC[_qp] * _T_grad[_qp](_component) -
          _sbC[_qp] * _T[_qp] * _ecC[_qp] * _potential_E_int_grad[_qp](_component) -
          _sbC[_qp] * _sbC[_qp] * _ecC[_qp] * _T[_qp] * _T_grad[_qp](_component)) *
         _len_scale;
}
