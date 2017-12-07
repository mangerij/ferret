#ifndef POLARIZATIONCOMPONENTVALUE_H
#define POLARIZATIONCOMPONENTVALUE_H

#include "ElementIntegralPostprocessor.h"

class PolarizationComponentValue;

template<>
InputParameters validParams<PolarizationComponentValue>();

class PolarizationComponentValue: public ElementIntegralPostprocessor
{
public:
  PolarizationComponentValue(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const VariableValue & _polar;
};

#endif
