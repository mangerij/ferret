/**
 * @file   TotalEnergySkFlow.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Aug 15 15:54:15 2013
 *
 * @brief This is a total energy postprocessor that tracks the following
 *        bulk, wall, coupled elastic, depolarization, and electrostatic energy. 
 *        The term "flow" means the energy that is nonzero under variational 
 *        differentiation in the diffusive gradient flow approach.
 */

#include "TotalEnergySkFlow.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergySkFlow>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0, "name of electrostatic energy postprocessor");
  params.addParam<PostprocessorName>("Fcoupled", 0.0, "name of the coupled energy postprocessor");
  params.addParam<PostprocessorName>("Fdepol", 0.0, "name of the depolarization energy postprocessor");
  return params;
}

TotalEnergySkFlow::TotalEnergySkFlow(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec"))),
  _Fcoupled(getPostprocessorValue(getParam<PostprocessorName>("Fcoupled"))),
  _Fdepol(getPostprocessorValue(getParam<PostprocessorName>("Fdepol")))
{
}

TotalEnergySkFlow::~TotalEnergySkFlow(){
}

void
TotalEnergySkFlow::initialize(){
}

void
TotalEnergySkFlow::execute(){
}

Real
TotalEnergySkFlow::getValue()
{
  return _Fbulk + _Fwall + _Felec + _Fcoupled + _Fdepol;
}
