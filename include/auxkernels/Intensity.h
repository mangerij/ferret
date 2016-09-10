#ifndef INTENSITY_H
#define INTENSITY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class Intensity;

template<>
InputParameters validParams<Intensity>();


class Intensity : public AuxKernel
{
public:
  Intensity(const InputParameters & parameters);

  virtual ~Intensity() {}

protected:
  virtual Real computeValue();
private:
  const VariableValue & _ReE_x;
  const VariableValue & _ReE_y;
  const VariableValue & _ReE_z;
  const VariableValue & _ImagE_x;
  const VariableValue & _ImagE_y;
  const VariableValue & _ImagE_z;
  //const VariableValue & _ReH_x;
  //const VariableValue & _ReH_y;
  //const VariableValue & _ReH_z;
  //const VariableValue & _ImagH_x;
 // const VariableValue & _ImagH_y;
 // const VariableValue & _ImagH_z;

};

#endif // INTENSITY_H
