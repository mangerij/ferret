/**
 * @file   TotalEnergyP.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Mar 3 2017
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGYP_H
#define TOTALENERGYP_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyP;

template<>
InputParameters validParams<TotalEnergyP>();

//TODO: change the base class!
class TotalEnergyP : public GeneralPostprocessor
{
public:
  TotalEnergyP(const InputParameters & parameters);
  virtual ~TotalEnergyP();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall, & _Fec, & _Fdepol;
};

#endif
