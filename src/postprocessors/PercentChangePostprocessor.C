
#include "PercentChangePostprocessor.h"

template<>
InputParameters validParams<PercentChangePostprocessor>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<PostprocessorName>("postprocessor", "The name of the postprocessor used for exit criterion");
  return params;
}

PercentChangePostprocessor::PercentChangePostprocessor(const InputParameters & parameters) :
    GeneralPostprocessor(parameters),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _postprocessor_old(getPostprocessorValueOld("postprocessor"))

{
}

void
PercentChangePostprocessor::initialize(){
}

void
PercentChangePostprocessor::execute(){
}

Real
PercentChangePostprocessor::getValue()
{
  return fabs( ( fabs(_postprocessor)- fabs(_postprocessor_old) )*pow(fabs(_postprocessor),-1));
}

