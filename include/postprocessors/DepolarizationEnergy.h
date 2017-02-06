/**
 * @file   DepolarizationEnergy.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief compute the electric depolarization energy
 *
 *
 */

#ifndef DEPOLARIZATIONENERGY_H
#define DEPOLARIZATIONENERGY_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class DepolarizationEnergy;

template<>
InputParameters validParams<DepolarizationEnergy>();

//TODO: change the base class!
class DepolarizationEnergy : public ElementIntegralPostprocessor
{
public:
  DepolarizationEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue & _polar_z;
  const Real _len_scale;
  const Real _avePz;
  const Real _lambda;
  const Real _permitivitty;
};

#endif
