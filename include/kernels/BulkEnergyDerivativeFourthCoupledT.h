/**
 * @file   BulkEnergyDerivativeFourthCoupledT.h
 * @author J. Mangeri <mangerij@anl.gov>
 * @date   Thu Aug 13 2:00 2015
 *
 */

#ifndef BULKENERGYDERIVATIVEFOURTHCOUPLEDT_H
#define BULKENERGYDERIVATIVEFOURTHCOUPLEDT_H

#include "Kernel.h"

class BulkEnergyDerivativeFourthCoupledT;

template<>
InputParameters validParams<BulkEnergyDerivativeFourthCoupledT>();

class BulkEnergyDerivativeFourthCoupledT: public Kernel
{
public:

  BulkEnergyDerivativeFourthCoupledT(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const unsigned int _temperature_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableValue & _temperature;
  const Real _alpha0, _alpha11, _alpha12, _Tc;
  const Real _len_scale;
};
#endif //BULKENERGYDERIVATIVEFOURTHCOUPLEDT_H
