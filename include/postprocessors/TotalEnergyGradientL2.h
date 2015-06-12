/**
 * @file   TotalEnergyGradientL2.h
 * @date   Thu Aug 15 15:48:51 2013
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGYGRADIENTL2_H
#define TOTALENERGYGRADIENTL2_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyGradientL2;

template<>
InputParameters validParams<TotalEnergyGradientL2>();

//TODO: change the base class!
class TotalEnergyGradientL2 : public GeneralPostprocessor
{
public:
  TotalEnergyGradientL2(const std::string & name, InputParameters parameters);
  virtual ~TotalEnergyGradientL2();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue &_gradx,&_grady,&_gradz;
};

#endif
