/**
 * @file   AdhocConstIC.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Aug 14 11:41:18 2013
 *
 * @brief
 *
 *
 */
#include "AdhocConstIC.h"

#include "libmesh/point.h"
#include <limits>
#include <cmath>
template<>
InputParameters validParams<AdhocConstIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("value0", "value for z<0.5");
  params.addRequiredParam<Real>("value1", "value for z>=0.5");
  return params;
}



AdhocConstIC::AdhocConstIC(const InputParameters & parameters) :
  InitialCondition(parameters),
  _val0(getParam<Real>("value0")),
  _val1(getParam<Real>("value1"))
{
}

Real
AdhocConstIC::value(const Point & p)
{
  if(p(2)-0.5<1000*std::numeric_limits<Real>::epsilon()) return _val0;
  else return _val1;
}
