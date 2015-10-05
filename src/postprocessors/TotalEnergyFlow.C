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
  params.addParam<PostprocessorName>("bulk_energy", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("wall_energy", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("bulk_energy_fourth", 0.0,  "name of bulk energy postprocessor to fourth order");
  params.addParam<PostprocessorName>("electrostatic_energy", 0.0, "name of electrostatic energy postprocessor");
  params.addParam<PostprocessorName>("coupled_energy", 0.0, "name of the coupled energy postprocessor");
  return params;
}

TotalEnergyFlow::TotalEnergyFlow(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _bulk_energy(getPostprocessorValue(getParam<PostprocessorName>("bulk_energy"))),
  _wall_energy(getPostprocessorValue(getParam<PostprocessorName>("wall_energy"))),
  _bulk_energy_fourth(getPostprocessorValue(getParam<PostprocessorName>("bulk_energy_fourth"))),
  _electrostatic_energy(getPostprocessorValue(getParam<PostprocessorName>("electrostatic_energy"))),
  _coupled_energy(getPostprocessorValue(getParam<PostprocessorName>("coupled_energy")))
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
  return _bulk_energy + _wall_energy + _bulk_energy_fourth + _electrostatic_energy + _coupled_energy;
}
