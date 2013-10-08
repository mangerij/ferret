/**
 * @file   ElectrostaticEnergyDensityTotal.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 17:28:22 2013
 * @brief calculate \permittivity*|\grad\phi|^2+P\cdot\grad\phi
 */

#ifndef ELECTROSTATICENERGYDENSITYTOTAL_H
#define ELECTROSTATICENERGYDENSITYTOTAL_H

#include "AuxKernel.h"


//Forward Declarations
class ElectrostaticEnergyDensityTotal;

template<>
InputParameters validParams<ElectrostaticEnergyDensityTotal>();

/**
 * Coupled auxiliary value
 */
class ElectrostaticEnergyDensityTotal : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  ElectrostaticEnergyDensityTotal(const std::string & name, InputParameters parameters);

protected:
  const Real _permittivity;
  const VariableGradient& _potential_grad;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  virtual Real computeValue();
  const Real _len_scale;

};

#endif // ELECTROSTATICENERGYDENSITYTOTAL_H
