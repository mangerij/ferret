/**
 * @file   TotalEnergySkFlow.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGYSKFLOW_H
#define TOTALENERGYSKFLOW_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergySkFlow;

template<>
InputParameters validParams<TotalEnergySkFlow>();

//TODO: change the base class!
class TotalEnergySkFlow : public GeneralPostprocessor
{
public:
  TotalEnergySkFlow(const InputParameters & parameters);
  virtual ~TotalEnergySkFlow();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall, & _Felec, & _Fcoupled, & _Fdepol;
};

#endif
