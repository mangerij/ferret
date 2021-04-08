#include "q1.h"
#include "Material.h"

registerMooseObject("FerretApp", q1);

template <>
InputParameters
validParams<q1>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to k*deltaT = 0");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction the variable "
                                        "this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<MaterialPropertyName>(
      "thC", "Thermal Conductivity", "Property name of the thermal conductivity material property");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

q1::q1(const InputParameters & parameters)
  : Kernel(parameters),
    _thC(getMaterialProperty<Real>("thC")),
    _component(getParam<unsigned int>("component")),
    _T_var(coupled("T")),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _len_scale(getParam<Real>("len_scale"))
{
}

Real
q1::computeQpResidual()
{
  return -_grad_test[_i][_qp](_component) * (_thC[_qp] * _T_grad[_qp](_component)) * _len_scale;
}

Real
q1::computeQpJacobian()
{
  return -_grad_test[_i][_qp](_component) * (_thC[_qp] * _grad_phi[_j][_qp](_component)) *
         _len_scale;
}
