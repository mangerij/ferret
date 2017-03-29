/**
 * @file   TotalEnergyG.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Mar 3 2017
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGYG_H
#define TOTALENERGYG_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyG;

template<>
InputParameters validParams<TotalEnergyG>();

//TODO: change the base class!
class TotalEnergyG : public GeneralPostprocessor
{
public:
  TotalEnergyG(const InputParameters & parameters);
  virtual ~TotalEnergyG();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall, & _Faniso, & _Fec;
};

#endif
