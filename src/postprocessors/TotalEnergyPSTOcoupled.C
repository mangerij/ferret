

#include "TotalEnergyPSTOcoupled.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyPSTOcoupled>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("FbulkPSTO", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0,  "name of electrical energy postprocessor");
  params.addParam<PostprocessorName>("Fcoupled", 0.0,  "name of coupled energy postprocessor");

  return params;
}

TotalEnergyPSTOcoupled::TotalEnergyPSTOcoupled(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _FbulkPSTO(getPostprocessorValue(getParam<PostprocessorName>("FbulkPSTO"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec"))),
  _Fcoupled(getPostprocessorValue(getParam<PostprocessorName>("Fcoupled")))


{
}

TotalEnergyPSTOcoupled::~TotalEnergyPSTOcoupled(){
}

void
TotalEnergyPSTOcoupled::initialize(){
}

void
TotalEnergyPSTOcoupled::execute(){
}

Real
TotalEnergyPSTOcoupled::getValue()
{
  ///  return _bulk_energy + _wall_energy + _electrical_energy + _coupled_energy +;
  return _FbulkPSTO + _Fwall + _Felec + _Fcoupled;
}
