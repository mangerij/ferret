/**
 * @file   PiezoelectricApprox.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate piezoelectric constant 
 * d_{33} = P_z / \sigma_{zz}
 * for more information, see IEEE Trans. Ultra. Ferro. Freq. Contr. 51, 3, (2004)
 *
 */


#ifndef PIEZOELECTRICAPPROX_H
#define PIEZOELECTRICAPPROX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class PiezoelectricApprox;

template<>
InputParameters validParams<PiezoelectricApprox>();


class PiezoelectricApprox : public AuxKernel
{
public:
  PiezoelectricApprox(const InputParameters & parameters);

  virtual ~PiezoelectricApprox() {}

protected:
  virtual Real computeValue();

private:
  const MaterialProperty<RankTwoTensor> & _stress;
  const VariableValue & _polar_z;
};

#endif // PIEZOELECTRICAPPROX_H
