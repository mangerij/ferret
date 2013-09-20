/**
 * @file   TotalEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:54:15 2013
 *
 * @brief
 *
 *
 */

#include "TotalEnergy.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergy>()
{
  //TODO: inherit from an appropriate postprocessor
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<PostprocessorName>("bulk_energy","name of bulk_energy postprocessor");
  params.addRequiredParam<PostprocessorName>("wall_energy","name of wall_energy postprocessor");
  params.addRequiredParam<PostprocessorName>("electric_energy","name of electric_energy postprocessor");
  return params;
}

TotalEnergy::TotalEnergy(const std::string & name, InputParameters parameters) :
  GeneralPostprocessor(name, parameters),
  _bulk_energy(getPostprocessorValue(getParam<PostprocessorName>("bulk_energy"))),
  _wall_energy(getPostprocessorValue(getParam<PostprocessorName>("wall_energy"))),
  _electric_energy(getPostprocessorValue(getParam<PostprocessorName>("electric_energy")))
{
}

TotalEnergy::~TotalEnergy(){
}

void
TotalEnergy::initialize(){
}

void
TotalEnergy::execute(){
}

Real
TotalEnergy::getValue()
{
  return _bulk_energy+_wall_energy+_electric_energy;
}
