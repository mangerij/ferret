/**
 * @file   BulkEnergyDerivativePSTO.h
 * @author J. Mangeri <john.mangeri@uconn.edu> and S. Churchill <steven.churchill@uconn.edu>
 * @date   Jun 1 12:00:20 2015
 *
 * Bulk energy function derived from first-principles calculations 
 * for more information, see Mangeri et al npj Computational Materials 2, 16020 (2016)
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

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;

  const Real _alpha1;
  const Real _alpha2;
  const Real _alpha3;
  const Real _alpha4;
  const Real _alpha5;
  const Real _x1;
  const Real _x2;
  const Real _x3;
  const Real _x4;
  const Real _x5;
  const Real _x6;
  const Real _epsilon;
  const Real _T;
};
#endif //BULKENERGYDERIVATIVEPSTO_H
