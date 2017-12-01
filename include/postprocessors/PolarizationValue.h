#ifndef POLARIZATIONVALUE_H
#define POLARIZATIONVALUE_H

#include "ElementIntegralPostprocessor.h"

class PolarizationValue;

template<>
InputParameters validParams<PolarizationValue>();

class PolarizationValue: public ElementIntegralPostprocessor
{
public:
  PolarizationValue(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
};

#endif
