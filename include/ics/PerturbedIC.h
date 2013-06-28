/**
 * @file   PerturbedIC.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Jun 26 18:40:34 2013
 *
 * @brief
 *
 *
 */
#ifndef PERTURBEDIC_H
#define PERTURBEDIC_H

#include "Kernel.h"
#include "InitialCondition.h"
#include "InputParameters.h"

// Forward Declarations
class PerturbedIC;

namespace libMesh { class Point; }

template<>
InputParameters validParams<PerturbedIC>();


class PerturbedIC:public InitialCondition
{
public:
  PerturbedIC(const std::string & name, InputParameters parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);
private:
  Real _mean;
  Real _factor;
};

#endif //PerturbedIC
