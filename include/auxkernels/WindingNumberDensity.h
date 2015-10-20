/**
 * @file   WindingNumberDensity.h
 * @brief Calculates the topological winding number density in 3D for
 *        the polarization vector P.
 */

#ifndef WINDINGNUMBERDENSITY_H
#define WINDINGNUMBERDENSITY_H

#include "AuxKernel.h"


//Forward Declarations
class WindingNumberDensity;

template<>
InputParameters validParams<WindingNumberDensity>();

class WindingNumberDensity : public AuxKernel
{
public:
  WindingNumberDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
};

#endif // WINDINGNUMBERDENSITY_H
