/**
 * @file   TotalEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Aug 15 15:48:51 2013
 *
 * @brief
 *
 *
 */

#ifndef TOTALENERGY_H
#define TOTALENERGY_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergy;

template<>
InputParameters validParams<TotalEnergy>();

//TODO: change the base class!
class TotalEnergy : public GeneralPostprocessor
{
public:
  TotalEnergy(const InputParameters & parameters);
  virtual ~TotalEnergy();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _bulk_energy, & _wall_energy,& _bulk_energy_fourth, & _electrostatic_energy, & _elastic_energy, & _coupled_energy;
};

#endif
