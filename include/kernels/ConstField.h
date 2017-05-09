/**
 * @file   ConstField.h
 * @author J. Mangeri <john.mangeri@uconn.edu
 */

#ifndef CONSTFIELD_H
#define CONSTFIELD_H

#include "Kernel.h"

//Forward Declarations
class ConstField;

template<>
InputParameters validParams<ConstField>();

class ConstField: public Kernel
{
public:

  ConstField(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned int _polar_z_var;
  const VariableValue & _polar_z;

  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const Real _field;

};
#endif //CONSTFIELD_H
