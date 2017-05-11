

#include "TotalEnergyPSTO.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyPSTO>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("FbulkPSTO", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0,  "name of electrical energy postprocessor");

  return params;
}

TotalEnergyPSTO::TotalEnergyPSTO(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _FbulkPSTO(getPostprocessorValue(getParam<PostprocessorName>("FbulkPSTO"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec")))


{
}

TotalEnergyPSTO::~TotalEnergyPSTO(){
}

void
TotalEnergyPSTO::initialize(){
}

void
TotalEnergyPSTO::execute(){
}

Real
TotalEnergyPSTO::getValue()
{
  ///  return _bulk_energy + _wall_energy + _electrical_energy;
  return _FbulkPSTO + _Fwall + _Felec;
}
