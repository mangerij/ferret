/**
 * @file   TotalEnergyFlowNoElast.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Aug 15 15:54:15 2013
 *
 * @brief
 *
 *
 */

#include "TotalEnergyFlowNoElast.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyFlowNoElast>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0, "name of electrostatic energy postprocessor");
  return params;
}

TotalEnergyFlowNoElast::TotalEnergyFlowNoElast(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec")))
{
}

TotalEnergyFlowNoElast::~TotalEnergyFlowNoElast(){
}

void
TotalEnergyFlowNoElast::initialize(){
}

void
TotalEnergyFlowNoElast::execute(){
}

Real
TotalEnergyFlowNoElast::getValue()
{
  ///  return _bulk_energy + _wall_energy + _electrostatic_energy;
  return _Fbulk + _Fwall + _Felec;
}
