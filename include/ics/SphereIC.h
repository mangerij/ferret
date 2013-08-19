/**
 * @file   SphereIC.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Aug 12 16:46:06 2013
 *
 * @brief
 *
 *
 */
#ifndef SPHEREIC_H
#define SPHEREIC_H

#include "Kernel.h"
#include "InitialCondition.h"
#include "InputParameters.h"
// Forward Declarations
class SphereIC;
class Function;

namespace libMesh { class Point; }

template<>
InputParameters validParams<SphereIC>();


class SphereIC:public InitialCondition
{
public:
  SphereIC(const std::string & name, InputParameters parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);
protected:
  Function & _radial_func;
  Function & _polar_func;
  Function & _azimuthal_func;
  unsigned int _index;
};

#endif //SphereIC
