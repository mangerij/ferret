#include "SeebeckEffect.h"

registerMooseObject("FerretApp", SeebeckEffect);

template <>
InputParameters
validParams<SeebeckEffect>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a contribution due to nabla.j = 0");
  params.addParam<MaterialPropertyName>("sbC", "seebeck_coefficient");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction the variable "
                                        "this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("potential_E_int", "electrical potential");
  params.addRequiredCoupledVar("T", "temperature");
  // params.addParam<MaterialPropertyName>("ecC", "Electrical Conductivity", "Property name of the
  // electrical conductivity material property");
  params.addParam<MaterialPropertyName>(
      "sbC", "Seebeck coefficient", "Property name of the Seebeck coefficient material property");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

SeebeckEffect::SeebeckEffect(const InputParameters & parameters)
  : Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _potential_E_int_var(coupled("potential_E_int")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _T_var(coupled("T")),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    // _ecC(getMaterialProperty<Real>("ecC")),
    _sbC(getMaterialProperty<Real>("sbC")),
    _len_scale(getParam<Real>("len_scale"))
{
}

Real
SeebeckEffect::computeQpResidual()
{
  return -_grad_test[_i][_qp](_component) *
         (_potential_E_int_grad[_qp](_component) + _sbC[_qp] * _T_grad[_qp](_component)) *
         _len_scale;
}

Real
SeebeckEffect::computeQpJacobian()
{
  // return - _grad_test[_i][_qp](_component) * (_potential_E_int_grad[_qp](_component) + _sbC[_qp]
  // * _grad_phi[_j][_qp](_component)) * _len_scale;
  return -_grad_test[_i][_qp](_component) *
         (_grad_phi[_j][_qp](_component) + _sbC[_qp] * _T_grad[_qp](_component)) * _len_scale;
}

Real
SeebeckEffect::computeQpOffDiagJacobian(unsigned int jvar)
{
  // if(jvar == _potential_E_int_var)
  //  {
  //   return (-_grad_test[_i][_qp](_component)) * _grad_phi[_j][_qp](_component) * _len_scale;
  //
  //   }
  //   else
  //   {
  //     return 0.0;
  //   }
  if (jvar == _T_var)
  {
    return (-_grad_test[_i][_qp](_component)) * _sbC[_qp] * _grad_phi[_j][_qp](_component) *
           _len_scale;
  }
  else
  {
    return 0.0;
  }
}
