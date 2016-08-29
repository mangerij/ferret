#include "AngleAux.h"

template<>

InputParameters validParams<AngleAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}


AngleAux::AngleAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{
}

Real
AngleAux::computeValue()
{
//(180 / libMesh::pi )
    return  std::acos(_polar_z[_qp] / std::pow(std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0) + std::pow(_polar_z[_qp], 2.0), 0.5));
}


