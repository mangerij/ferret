#ifndef SCREENAUX_H
#define SCREENAUX_H

#include "AuxKernel.h"


//Forward declarations
class ScreenAux;

template<>
InputParameters validParams<ScreenAux>();

class ScreenAux : public AuxKernel
{
public:
  ScreenAux( const std::string & name, InputParameters parameters );

  virtual ~ScreenAux() {}

protected:
  virtual Real computeValue();
  const MooseArray<Point> & _normals;

private:
  const VariableGradient &  _potential_int_grad;

};

#endif
