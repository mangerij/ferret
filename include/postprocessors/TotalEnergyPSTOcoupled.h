

#ifndef TOTALENERGYPSTOCOUPLED_H
#define TOTALENERGYPSTOCOUPLED_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyPSTOcoupled;

template<>
InputParameters validParams<TotalEnergyPSTOcoupled>();

//TODO: change the base class!
class TotalEnergyPSTOcoupled : public GeneralPostprocessor
{
public:
  TotalEnergyPSTOcoupled(const InputParameters & parameters);
  virtual ~TotalEnergyPSTOcoupled();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _FbulkPSTO, & _Fwall, & _Felec, & _Fcoupled;
};

#endif
