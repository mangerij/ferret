#ifndef CURLPMAG_H
#define CURLPMAG_H

#include "AuxKernel.h"

//Forward declarations
class CurlPMag;

template<>
InputParameters validParams<CurlPMag>();

class CurlPMag : public AuxKernel
{
public:
  CurlPMag(const InputParameters & parameters);

  virtual ~CurlPMag() {}

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
