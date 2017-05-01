#include "PkNorm.h"
#include <math.h>

template<>
InputParameters validParams<PkNorm>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("component", "the component of the normalized vector to store");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}


PkNorm::PkNorm(const InputParameters & parameters) :
  AuxKernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{}

Real
PkNorm::computeValue()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  if (_component == 0)
    return _polar_x[_qp] / sqrt(w*w);
  else if (_component == 1)
    return _polar_y[_qp] / sqrt(w*w);
  else if (_component == 2)
    return _polar_z[_qp] / sqrt(w*w);
  else
    return 0.0;
}
