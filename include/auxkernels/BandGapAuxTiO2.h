#ifndef BANDGAPAUXTIO2_H
#define BANDGAPAUXTIO2_H

#include "AuxKernel.h"
#include "TensorMechanicsMaterial.h" //may not need this
#include "RankTwoTensor.h"

//Forward declarations
class BandGapAuxTiO2;

template<>
InputParameters validParams<BandGapAuxTiO2>();


class BandGapAuxTiO2 : public AuxKernel
{
public:
  BandGapAuxTiO2(const InputParameters & parameters);

  virtual ~BandGapAuxTiO2() {}

protected:
  virtual Real computeValue();

private:
  const MaterialProperty<RankTwoTensor> & _stress;
  const Real _ba, _bc, _E0;

};

#endif // BANDGAP_H
