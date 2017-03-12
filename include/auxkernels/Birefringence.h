/**
 * @file   Birefringence.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate the birefringence
 * \Delta n = n_o - n_e
 *
 */

#ifndef BIREFRINGENCE_H
#define BIREFRINGENCE_H

#include "AuxKernel.h"

//Forward declarations
class Birefringence;

template<>
InputParameters validParams<Birefringence>();


class Birefringence : public AuxKernel
{
public:
  Birefringence(const InputParameters & parameters);

  virtual ~Birefringence() {}

protected:
  virtual Real computeValue();
  const VariableValue & _var1;
  const VariableValue & _var2;

private:

};

#endif // BIREFRINGENCE_H
