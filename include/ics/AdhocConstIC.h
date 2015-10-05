/**
 * @file   AdhocConstIC.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Aug 14 11:39:24 2013
 *
 * @brief provide piecewise constant: if z<0.5, one value otherwise another value
 *
 *
 */
#ifndef ADHOCCONSTIC_H
#define ADHOCCONSTIC_H

#include "Kernel.h"
#include "InitialCondition.h"
#include "InputParameters.h"

// Forward Declarations
class AdhocConstIC;

namespace libMesh { class Point; }

template<>
InputParameters validParams<AdhocConstIC>();


class AdhocConstIC:public InitialCondition
{
public:
  AdhocConstIC(const InputParameters & parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);
protected:
  Real _val0,_val1;
};

#endif //AdhocConstIC
