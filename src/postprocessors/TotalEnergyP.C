/**
 * @file   TotalEnergyP.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Mar 3 2017
 *
 * @brief
 *
 *
 */

#include "TotalEnergyP.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyP>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Fec", 0.0,  "name of elec energy postprocessor");
  params.addParam<PostprocessorName>("Fdepol", 0.0,  "name of depol energy postprocessor");
  return params;
}

TotalEnergyP::TotalEnergyP(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Fec(getPostprocessorValue(getParam<PostprocessorName>("Fec"))),
  _Fdepol(getPostprocessorValue(getParam<PostprocessorName>("Fdepol")))
{
}

TotalEnergyP::~TotalEnergyP(){
}

void
TotalEnergyP::initialize(){
}

void
TotalEnergyP::execute(){
}

Real
TotalEnergyP::getValue()
{
  ///  return _bulk_energy + _wall_energy + _anisotropy_energy + _elec_energy;
  return _Fbulk + _Fwall + _Fec + _Fdepol;
}
