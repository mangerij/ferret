/**
 * @file   DepolEnergy.h
 * @author J. Mangeri <john.mangeri@uconn.edu
 */

#ifndef DEPOLENERGY_H
#define DEPOLENERGY_H

#include "Kernel.h"

//Forward Declarations
class DepolEnergy;

template<>
InputParameters validParams<DepolEnergy>();

class DepolEnergy: public Kernel
{
public:

  DepolEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned int _polar_z_var;
  const VariableValue & _polar_z;

  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const PostprocessorValue & _avePz;
  const Real _lambda;
  const Real _permitivitty;

};
#endif //DEPOLENERGY_H
