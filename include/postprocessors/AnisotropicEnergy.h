/**
 * @file   AnisotropicEnergy.h
 * @author J. Mangeri <john.mangeri@uconn.edu
 */

#ifndef ANISOTROPICENERGY_H
#define ANISOTROPICENERGY_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class AnisotropicEnergy;

template<>
InputParameters validParams<AnisotropicEnergy>();

class AnisotropicEnergy: public ElementIntegralPostprocessor
{
public:

  AnisotropicEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

private:
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const Real _K;

};
#endif //ANISOTROPICENERGY_H
