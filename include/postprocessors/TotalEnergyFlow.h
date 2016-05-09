/**
 * @file   TotalEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:48:51 2013
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGYFLOW_H
#define TOTALENERGYFLOW_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyFlow;

template<>
InputParameters validParams<TotalEnergyFlow>();

//TODO: change the base class!
class TotalEnergyFlow : public GeneralPostprocessor
{
public:
  TotalEnergyFlow(const InputParameters & parameters);
  virtual ~TotalEnergyFlow();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall, & _Felec, & _Fcoupled;
};

#endif
