/**
 * @file   TotalEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:54:15 2013
 *
 * @brief
 *
 *
 */

#include "TotalEnergyGradientL2.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyGradientL2>()
{
  //TODO: inherit from an appropriate postprocessor
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<PostprocessorName>("gradx","name of energy x gradient postprocessor");
  params.addRequiredParam<PostprocessorName>("grady","name of energy y gradient postprocessor");
  params.addRequiredParam<PostprocessorName>("gradz","name of energy z gradient postprocessor");
  return params;
}

TotalEnergyGradientL2::TotalEnergyGradientL2(const std::string & name, InputParameters parameters) :
  GeneralPostprocessor(name, parameters),
  _gradx(getPostprocessorValue(getParam<PostprocessorName>("gradx"))),
  _grady(getPostprocessorValue(getParam<PostprocessorName>("grady"))),
  _gradz(getPostprocessorValue(getParam<PostprocessorName>("gradz")))
{
}

TotalEnergyGradientL2::~TotalEnergyGradientL2(){
}

void
TotalEnergyGradientL2::initialize(){
}

void
TotalEnergyGradientL2::execute(){
}

Real
TotalEnergyGradientL2::getValue()
{
  return pow(pow(_gradx,2)+pow(_grady,2)+pow(_gradz,2),0.5);
}
