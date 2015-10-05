#include "TensorPressureAux.h"

template<>
InputParameters validParams<TensorPressureAux>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

TensorPressureAux::TensorPressureAux(const InputParameters & parameters) :
  AuxKernel(parameters ),
   _stress( getMaterialProperty<RankTwoTensor>("stress") )
{}

Real
TensorPressureAux::computeValue()
{
    return -0.33333333333*_stress[_qp].trace();
    
}


