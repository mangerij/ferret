/**
 * @file   BulkEnergyDerivativePSTO.h
 * @author J. Mangeri <john.mangeri@uconn.edu> and S. Churchill <steven.churchill@uconn.edu>
 * @date   Jun 1 12:00:20 2015
 *
 * @brief
 *
 *
 */

#ifndef BULKENERGYDERIVATIVEPSTO_H
#define BULKENERGYDERIVATIVEPSTO_H

#include "Kernel.h"

class BulkEnergyDerivativePSTO;

template<>
InputParameters validParams<BulkEnergyDerivativePSTO>();

class BulkEnergyDerivativePSTO: public Kernel
{
public:

  BulkEnergyDerivativePSTO(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;

  const Real _alpha1, _alpha2, _alpha3, _alpha4, _alpha5,_x1, _x2, _x3, _x4, _x5, _x6, _epsilon, _T, _Tc;
};
#endif //BULKENERGYDERIVATIVEPSTO_H
