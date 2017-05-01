/**
 * @file   TotalWinding.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   
 *
 * @brief Integral over the volume of the winding number density q.
 *
 */


#ifndef TOTALWINDING_H
#define TOTALWINDING_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class TotalWinding;

template<>
InputParameters validParams<TotalWinding>();

//TODO: change the base class!
class TotalWinding : public ElementIntegralPostprocessor
{
public:
  TotalWinding(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue & _q;

};

#endif
