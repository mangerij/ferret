/**
 * @file   TotalEnergyFlowNoElastNoElec.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Aug 15 15:54:15 2013
 *
 * @brief
 *
 *
 */

#include "TotalEnergyFlowNoElastNoElec.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyFlowNoElastNoElec>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  return params;
}

TotalEnergyFlowNoElastNoElec::TotalEnergyFlowNoElastNoElec(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall")))
{
}

TotalEnergyFlowNoElastNoElec::~TotalEnergyFlowNoElastNoElec(){
}

void
TotalEnergyFlowNoElastNoElec::initialize(){
}

void
TotalEnergyFlowNoElastNoElec::execute(){
}

Real
TotalEnergyFlowNoElastNoElec::getValue()
{
  ///  return _bulk_energy + _wall_energy;
  return _Fbulk + _Fwall;
}
