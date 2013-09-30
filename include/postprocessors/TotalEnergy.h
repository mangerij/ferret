/**
 * @file   TotalEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:48:51 2013
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGY_H
#define TOTALENERGY_H

#include "GeneralPostprocessor.h"
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class TotalEnergy;

template<>
InputParameters validParams<TotalEnergy>();

//TODO: change the base class!
class TotalEnergy : public GeneralPostprocessor
{
public:
  TotalEnergy(const std::string & name, InputParameters parameters);
  virtual ~TotalEnergy();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue &_bulk_energy,&_wall_energy,&_electric_energy;
};

#endif
