/**
 * @file   TotalEnergyFlow.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Aug 15 15:54:15 2013
 *
 * @brief This is a total energy postprocessor that tracks the following
 *        bulk, wall, coupled elastic, and electrostatic energy. The term 
 *        "flow" means the energy that is nonzero under variational 
 *        differentiation in the diffusive gradient flow approach.
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
  return _Fbulk + _Fwall + _Felec + _Fcoupled;
}
