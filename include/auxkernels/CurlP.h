#ifndef CURLP_H
#define CURLP_H

#include "AuxKernel.h"

//Forward declarations
class CurlP;

template<>
InputParameters validParams<CurlP>();

class CurlP : public AuxKernel
{
public:
  CurlP(const InputParameters & parameters);

  virtual ~CurlP() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
};

#endif
