/**
 * @file   BulkEnergy.C
 * @author S. Churchell <steve.churchill@uconn.edu>
 *
 * @brief
 *  Total free energy postprocessor for PSTO (2D ferroelectric) 
 *
 */


#ifndef BULKENERGYPSTO_H
#define BULKENERGYPSTO_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class BulkEnergyPSTO;

template<>
InputParameters validParams<BulkEnergyPSTO>();

class BulkEnergyPSTO : public ElementIntegralPostprocessor
{
public:
  BulkEnergyPSTO(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const Real _alpha1, _alpha2, _alpha3, _alpha4, _alpha5,_x1, _x2, _x3, _x4, _x5, _x6, _epsilon, _T, _Tc;

};

#endif
