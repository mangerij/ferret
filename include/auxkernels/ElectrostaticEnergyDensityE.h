/**
 * @file   ElectrostaticEnergyDensityE.h   ElectrostaticEnergyDensityE.h
 * @author S. Gu <sgu@anl.gov>
 * @brief  Calculate ElectrostaticEnergyDensity |\grad\phi|^2
 */
/**
 * @file
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 14:51:01 2013
 *
 * @brief
 *
 *
 */

#ifndef ELECTROSTATICENERGYDENSITYE_H
#define ELECTROSTATICENERGYDENSITYE_H

#include "AuxKernel.h"


//Forward Declarations
class ElectrostaticEnergyDensityE;

template<>
InputParameters validParams<ElectrostaticEnergyDensityE>();

/**
 * Coupled auxiliary value
 */
class ElectrostaticEnergyDensityE : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  ElectrostaticEnergyDensityE(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  const VariableGradient & _potential_grad;
  //Real _value;
};

#endif // ELECTROSTATICENERGYDENSITYE_H
