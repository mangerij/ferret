/**
 * @file   TotalEnergyNoElast.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:48:51 2013
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGYFLOWNOELAST_H
#define TOTALENERGYFLOWNOELAST_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyFlowNoElast;

template<>
InputParameters validParams<TotalEnergyFlowNoElast>();

//TODO: change the base class!
class TotalEnergyFlowNoElast : public GeneralPostprocessor
{
public:
  TotalEnergyFlowNoElast(const InputParameters & parameters);
  virtual ~TotalEnergyFlowNoElast();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall, & _Felec;
};

#endif
