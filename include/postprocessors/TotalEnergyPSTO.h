

#ifndef TOTALENERGYPSTO_H
#define TOTALENERGYPSTO_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyPSTO;

template<>
InputParameters validParams<TotalEnergyPSTO>();

//TODO: change the base class!
class TotalEnergyPSTO : public GeneralPostprocessor
{
public:
  TotalEnergyPSTO(const InputParameters & parameters);
  virtual ~TotalEnergyPSTO();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _FbulkPSTO, & _Fwall, & _Felec;
};

#endif
