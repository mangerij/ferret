#include "CarrierInt.h"

registerMooseObject("FerretApp", CarrierInt);

template<>
InputParameters validParams<CarrierInt>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("ni", "intrinsic carriers");
  return params;
}

CarrierInt::CarrierInt(const InputParameters & parameters)
  :Kernel(parameters),
   _ni(getParam<Real>("ni"))
{
}

Real
CarrierInt::computeQpResidual()
{
  Real Rint = 0.0;
  Rint +=  _test[_i][_qp]* std::pow(_ni,2);

  return Rint;
}
Real
CarrierInt::computeQpJacobian()
{
  return 0.0;
}
