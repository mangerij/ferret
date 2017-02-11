/**
 * @file   SemiconductorChargeCarriers.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief
 *
 *
 */

#ifndef SEMICONDUCTORCHARGECARRIERS_H
#define SEMICONDUCTORCHARGECARRIERS_H

#include "Kernel.h"

class SemiconductorChargeCarriers;

template<>
InputParameters validParams<SemiconductorChargeCarriers>();

class SemiconductorChargeCarriers: public Kernel
{
public:

  SemiconductorChargeCarriers(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
   const unsigned int _potential_int_var;
   const VariableValue & _potential_int;
   const Real _q;
   const Real _kT;
   const Real _NA;
   const Real _NC;
   const Real _NV;
   const Real _EA;
   const Real _EC;
   const Real _EV;
   const Real _EF;
   const Real _len_scale;

};
#endif //SEMICONDUCTORCHARGECARRIERS_H
