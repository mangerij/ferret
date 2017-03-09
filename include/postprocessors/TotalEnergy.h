/**
 * @file   TotalEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Mar 3 2017
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGY_H
#define TOTALENERGY_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergy;

template<>
InputParameters validParams<TotalEnergy>();

//TODO: change the base class!
class TotalEnergy : public GeneralPostprocessor
{
public:
  TotalEnergy(const InputParameters & parameters);
  virtual ~TotalEnergy();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall;
};

#endif
