/**
 * @file   FluctuationKernel.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Jun 14 12:00:20 2016
 *
 * @brief This Kernel is used to introduce noise in between electric field changes in
 *        a quasi-static hysteresis calculation (see arxiv.org/pdf/1701.02613.pdf)
 */

#ifndef FLUCTUATIONKERNEL_H
#define FLUCTUATIONKERNEL_H

#include "Kernel.h"

class FluctuationKernel;

template<>
InputParameters validParams<FluctuationKernel>();

class FluctuationKernel: public Kernel
{
public:

  FluctuationKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const VariableValue & _deltaPi;
  const Real _len_scale;
};
#endif //FLUCTUATIONKERNEL_H
