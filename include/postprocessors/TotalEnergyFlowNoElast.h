/**
 * @file   TotalEnergyNoElast.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
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
