/**
 * @file   TotalEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:48:51 2013
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGYFLOW_H
#define TOTALENERGYFLOW_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyFlow;

template<>
InputParameters validParams<TotalEnergyFlow>();

//TODO: change the base class!
class TotalEnergyFlow : public GeneralPostprocessor
{
public:
  TotalEnergyFlow(const std::string & name, InputParameters parameters);
  virtual ~TotalEnergyFlow();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _bulk_energy, & _wall_energy,& _bulk_energy_fourth, & _electrostatic_energy, & _coupled_energy;
};

#endif
