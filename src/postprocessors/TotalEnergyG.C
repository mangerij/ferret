/**
 * @file   TotalEnergyG.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Mar 3 2017
 *
 * @brief
 *
 *
 */

#include "TotalEnergyG.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyG>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Faniso", 0.0,  "name of anisotropy energy postprocessor");
  params.addParam<PostprocessorName>("Fec", 0.0,  "name of elec energy postprocessor");
  return params;
}

TotalEnergyG::TotalEnergyG(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Faniso(getPostprocessorValue(getParam<PostprocessorName>("Faniso"))),
  _Fec(getPostprocessorValue(getParam<PostprocessorName>("Fec")))
{
}

TotalEnergyG::~TotalEnergyG(){
}

void
TotalEnergyG::initialize(){
}

void
TotalEnergyG::execute(){
}

Real
TotalEnergyG::getValue()
{
  ///  return _bulk_energy + _wall_energy + _anisotropy_energy + _elec_energy;
  return _Fbulk + _Fwall + _Faniso + _Fec;
}
