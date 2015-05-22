#ifndef BANDGAPAUXZNO_H
#define BANDGAPAUXZNO_H

#include "AuxKernel.h"
#include "TensorMechanicsMaterial.h" //may not need this
#include "RankTwoTensor.h"

//Forward declarations
class BandGapAuxZnO;

template<>
InputParameters validParams<BandGapAuxZnO>();


class BandGapAuxZnO : public AuxKernel
{
public:
  BandGapAuxZnO( const std::string & name, InputParameters parameters );

  virtual ~BandGapAuxZnO() {}

protected:
  virtual Real computeValue();

private:
  MaterialProperty<RankTwoTensor> & _strain;

private:

  Real _du;

private:

  Real _db;

private:

  Real _E0;

private:

  Real _Rb;

private:
  Real _nu;

};

#endif // BANDGAP_H


