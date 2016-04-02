/******************************************************
 *  Credit to A. Jokisaari
 *  Note: creates a flucutation about zero that is 
 *        spatially dependent but appears "random(ish)"
 *        and can be adjusted.
 *****************************************************/

#ifndef FLUCTUATIONSIC_H
#define FLUCTUATIONSIC_H

#include "InitialCondition.h"

// Forward Declarations
class FluctuationsIC;
namespace libMesh { class Point; }

template<>
InputParameters validParams<FluctuationsIC>();

class FluctuationsIC : public InitialCondition
{
public:
  FluctuationsIC(const InputParameters & parameters);
  virtual Real value(const Point & p);

protected:

private:
  Real _epsilon;
  Real _base_value;

  Point _q1;
  Point _q2;
  Point _q3;
  Point _q4;

  Real _h;

};

#endif //FLUCTUATIONSIC_H
