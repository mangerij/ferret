/**
 * @file   GrainSize.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Tue Jan 23 2017
 *
 * @brief
 *
 *
 */

#ifndef GRAINSIZE_H
#define GRAINSIZE_H


#include "GeneralPostprocessor.h"

//Forward Declarations
class GrainSize;

template<>
InputParameters validParams<GrainSize>();

class GrainSize : public GeneralPostprocessor
{
public:
  GrainSize(const InputParameters & parameters);
  virtual ~GrainSize();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _vol;
};

#endif
