/**
 * @file   TotalEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:54:15 2013
 *
 * @brief
 *
 *
 */

#include "TotalEnergyFlow.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyFlow>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0, "name of electrostatic energy postprocessor");
  params.addParam<PostprocessorName>("Fcoupled", 0.0, "name of the coupled energy postprocessor");
  return params;
}

TotalEnergyFlow::TotalEnergyFlow(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec"))),
  _Fcoupled(getPostprocessorValue(getParam<PostprocessorName>("Fcoupled")))
{
}

TotalEnergyFlow::~TotalEnergyFlow(){
}

void
TotalEnergyFlow::initialize(){
}

void
TotalEnergyFlow::execute(){
}

Real
TotalEnergyFlow::getValue()
{
  //  return _bulk_energy + _wall_energy + _electrostatic_energy;
  return _Fbulk + _Fwall + _Felec + _Fcoupled;
}
