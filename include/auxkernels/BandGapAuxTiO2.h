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
  BandGapAuxTiO2( const std::string & name, InputParameters parameters );

  virtual ~BandGapAuxTiO2() {}

protected:
  virtual Real computeValue();

private:
  MaterialProperty<RankTwoTensor> & _stress;

private:

  Real _ba;

private:

  Real _bc;

private:

  Real _E0;

};

#endif // BANDGAP_H


