/**
 * @file   TotalEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <john.mangeri@uconn.edu>
 * @date   Thu Mar 3 2017
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

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  return params;
}

TotalEnergy::TotalEnergy(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall")))
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
  ///  return _bulk_energy + _wall_energy;
  return _Fbulk + _Fwall;
}
