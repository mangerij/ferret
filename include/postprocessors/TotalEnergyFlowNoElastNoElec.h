/**
 * @file   TotalEnergyNoElastNoElec.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief
 *
 */

#ifndef TOTALENERGYFLOWNOELASTNOELEC_H
#define TOTALENERGYFLOWNOELASTNOELEC_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyFlowNoElastNoElec;

template<>
InputParameters validParams<TotalEnergyFlowNoElastNoElec>();

//TODO: change the base class!
class TotalEnergyFlowNoElastNoElec : public GeneralPostprocessor
{
public:
  TotalEnergyFlowNoElastNoElec(const InputParameters & parameters);
  virtual ~TotalEnergyFlowNoElastNoElec();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall;
};

#endif
