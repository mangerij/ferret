#include "Fourier_heat.h"
#include "Material.h"

registerMooseObject("FerretApp", Fourier_heat);

template<>
InputParameters validParams<Fourier_heat>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates a residual contribution due to k*deltaT = 0");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<MaterialPropertyName>("thC", "Thermal Conductivity", "Property name of the thermal conductivity material property");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

Fourier_heat::Fourier_heat(const InputParameters & parameters)
  :AuxKernel(parameters),
   _thC(getMaterialProperty<Real>("thC")),
   _component(getParam<unsigned int>("component")),
   _T_var(coupled("T")),
   _T(coupledValue("T")),
   _T_grad(coupledGradient("T")),
   _len_scale(getParam<Real>("len_scale"))
{
}
Real
Fourier_heat::computeValue()
{
  return -_thC[_qp] * _T_grad[_qp](_component) * _len_scale;
}
