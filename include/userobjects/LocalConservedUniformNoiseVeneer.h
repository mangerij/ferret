//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

/**
 * Veneer to build userobjects that generate a uniformly distributed random
 * number in the interval [-1:1] once per timestep for every quadrature point
 * in a way that the integral over all random numbers is zero.
 *
 * \see ConservedUniformNoise
 * \see ConservedMaskedUniformNoise
 */
template <class T>
class LocalConservedUniformNoiseVeneer : public T
{
public:
  LocalConservedUniformNoiseVeneer(const InputParameters & parameters);

protected:
  Real getQpRandom();
};

template <class T>
LocalConservedUniformNoiseVeneer<T>::LocalConservedUniformNoiseVeneer(const InputParameters & parameters)
  : T(parameters)
{
}

template <class T>
Real
LocalConservedUniformNoiseVeneer<T>::getQpRandom()
{
  return 2.0 * this->getRandomReal() - 1.0;
}
