
#ifndef PERCENTCHANGEPOSTPROCESSOR_H
#define PERCENTCHANGEPOSTPROCESSOR_H

#include "GeneralPostprocessor.h"

class PercentChangePostprocessor;

template<>
InputParameters validParams<PercentChangePostprocessor>();

/**
 * This postprocessor displays the change in the postprocessor between
 * adjacent timesteps
 */

 class PercentChangePostprocessor : public GeneralPostprocessor
 {
 public:
   PercentChangePostprocessor(const InputParameters & parameters);
   virtual void initialize();
   virtual void execute();
   virtual Real getValue();
 protected:
   const PostprocessorValue & _postprocessor, & _postprocessor_old;
 };

#endif /* PERCENTCHANGEPOSTPROCESSOR_H */
