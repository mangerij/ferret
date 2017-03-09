#ifndef SEMICONDUCINGCHARGECARRIERSAUX_H
#define SEMICONDUCINGCHARGECARRIERSAUX_H

#include "AuxKernel.h"


//Forward declarations
class SemiconductingChargeCarriersAux;

template<>
InputParameters validParams<SemiconductingChargeCarriersAux>();

class SemiconductingChargeCarriersAux : public AuxKernel
{
public:
  SemiconductingChargeCarriersAux(const InputParameters & parameters);

  virtual ~SemiconductingChargeCarriersAux() {}

protected:
  virtual Real computeValue();

private:
   const unsigned int _charge_type;
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
};

#endif /* SEMICONDUCINGCHARGECARRIERSAUX_H */
