/**
 * @file   Impermeability.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#ifndef IMPERMEABILITY_H
#define IMPERMEABILITY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class Impermeability;

template<>
InputParameters validParams<Impermeability>();


class Impermeability : public AuxKernel
{
public:
  Impermeability(const InputParameters & parameters);

  virtual ~Impermeability() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_i;
  const unsigned int _index_j;
  const MaterialProperty<RankTwoTensor> & _delta_beta_tensor;
  const MaterialProperty<RankTwoTensor> & _beta_tensor;

};

#endif // IMPERMEABILITY_H
