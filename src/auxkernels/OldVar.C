#include "OldVar.h"

template<>

InputParameters validParams<OldVar>()

{
  InputParameters params = validParams<AuxKernel>();
  return params;
}


OldVar::OldVar(const InputParameters & parameters) :
  AuxKernel(parameters)
{
}

Real
OldVar::computeValue()

{
    return _u_old[_qp];
}


