/**
 * @file   ThomasFermiTerm.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef THOMASFERMITERM_H
#define THOMASFERMITERM_H

#include "Kernel.h"

class ThomasFermiTerm;

template<>
InputParameters validParams<ThomasFermiTerm>();

class ThomasFermiTerm: public Kernel
{
public:

  ThomasFermiTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
   const unsigned int _potential_int_var;
   const VariableValue & _potential_int;
   const Real _len_scale;
   const Real _q;
   const Real _rho;
   const Real _EF;

};
#endif //THOMASFERMITERM_H
