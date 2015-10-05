/**
 * @file   BulkEnergyDerivativeSixth.h
 * @author J. Mangeri <mangerij@anl.gov>
 * @date   Thu Aug 13 2:00:20 2015
 *
 * @brief
 *
 *
 */

#ifndef BULKENERGYDERIVATIVESIXTHCOUPLEDT_H
#define BULKENERGYDERIVATIVESIXTHCOUPLEDT_H

#include "Kernel.h"

class BulkEnergyDerivativeSixthCoupledT;

template<>
InputParameters validParams<BulkEnergyDerivativeSixthCoupledT>();

class BulkEnergyDerivativeSixthCoupledT: public Kernel
{
public:

  BulkEnergyDerivativeSixthCoupledT(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const unsigned int _temperature_var;
  const VariableValue & _temperature;
  const Real _alpha0, _alpha11, _alpha12, _alpha111, _alpha112, _alpha123, _Tc;
  const Real _len_scale;
};
#endif //BULKENERGYDERIVATIVESIXTHCOUPLEDT_H
