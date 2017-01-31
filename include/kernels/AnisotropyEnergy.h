/**
 * @file   AnisotropyEnergy.h
 * @author J. Mangeri <john.mangeri@uconn.edu
 */

#ifndef ANISOTROPYENERGY_H
#define ANISOTROPYENERGY_H

#include "Kernel.h"

//Forward Declarations
class AnisotropyEnergy;

template<>
InputParameters validParams<AnisotropyEnergy>();

class AnisotropyEnergy: public Kernel
{
public:

  AnisotropyEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;

  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const Real _K;

};
#endif //ANISOTROPYENERGY_H
