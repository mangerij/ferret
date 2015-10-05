#ifndef TENSORPRESSUREAUX_H
#define TENSORPRESSUREAUX_H

#include "AuxKernel.h"
#include "TensorMechanicsMaterial.h" //may not need this
#include "RankTwoTensor.h"

//Forward declarations
class TensorPressureAux;

template<>
InputParameters validParams<TensorPressureAux>();


class TensorPressureAux : public AuxKernel
{
public:
  TensorPressureAux(const InputParameters & parameters);

  virtual ~TensorPressureAux() {}

protected:
  virtual Real computeValue();

private:
  const MaterialProperty<RankTwoTensor> & _stress;

};

#endif // TENSORPRESSUREAUX_H
