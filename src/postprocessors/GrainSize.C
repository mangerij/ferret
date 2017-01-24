/**
 * @file   GrainSize.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Tue Jan 23 2017
 *
 * @brief
 *
 *
 */

#include "GrainSize.h"
#include<iostream>

template<>
InputParameters validParams<GrainSize>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("vol", 0.0, "name of volume postprocessor");
  return params;
}

GrainSize::GrainSize(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _vol(getPostprocessorValue(getParam<PostprocessorName>("vol")))
{
}

GrainSize::~GrainSize(){
}

void
GrainSize::initialize(){
}

void
GrainSize::execute(){
}

Real
GrainSize::getValue()
{
  return 2.0 * std::pow(3.0 * _vol / (4.0 * 3.14159), 0.33333);
}
