/**
 * @file   WindingNumberDensity.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief calculate the winding number density in 3D for the polarization field P with 
 *        out-of-plane direction along z (this is the "t" coordinate system in 
 *        V. Stepkova and J. Hlinka, Phase Transitions, 2017 Vol. 90, No. 1, 11-16)
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
  const VariableValue & _norm_polar_x;
  const VariableValue & _norm_polar_y;
  const VariableValue & _norm_polar_z;
  const VariableGradient & _norm_polar_x_grad;
  const VariableGradient & _norm_polar_y_grad;
  const VariableGradient & _norm_polar_z_grad;
};

#endif // WINDINGNUMBERDENSITY_H
