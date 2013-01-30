#include "Magnetizing.h"

template<>
InputParameters validParams<Magnetizing>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

Magnetizing::Magnetizing(const std::string & name, InputParameters parameters) :
  Kernel(name, parameters),
  _p(getMaterialProperty<RealVectorValue>("Polarization"))
{
}

Real
Magnetizing::computeQpResidual()
{
   return -4.0*pi*_grad_test[_i][_qp]*_p[_qp];

}

