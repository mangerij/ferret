#ifndef SURFACECHARGEAUX_H
#define SURFACECHARGEAUX_H

#include "AuxKernel.h"


//Forward declarations
class SurfaceChargeAux;

template<>
InputParameters validParams<SurfaceChargeAux>();

class SurfaceChargeAux : public AuxKernel
{
public:
  SurfaceChargeAux( const std::string & name, InputParameters parameters );

  virtual ~SurfaceChargeAux() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
};

#endif 


