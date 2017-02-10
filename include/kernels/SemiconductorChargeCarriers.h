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

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
   const unsigned int _potential_int_var;
   const VariableValue & _potential_int;
   const Real _len_scale;

};
#endif //SEMICONDUCTORCHARGECARRIERS_H
