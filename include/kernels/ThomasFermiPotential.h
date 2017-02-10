/**
 * @file   ThomasFermiPotential.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef THOMASFERMIPOTENTIAL_H
#define THOMASFERMIPOTENTIAL_H

#include "Kernel.h"

class ThomasFermiPotential;

template<>
InputParameters validParams<ThomasFermiPotential>();

class ThomasFermiPotential: public Kernel
{
public:

  ThomasFermiPotential(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
   const unsigned int _potential_int_var;
   const VariableValue & _potential_int;
   const Real _len_scale;
   const Real _TFconstant;

};
#endif //THOMASFERMIPOTENTIAL_H
