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
  const MaterialProperty<RankTwoTensor> & _strain;
  const Real _du, _db, _E0, _Rb, _nu;

};

#endif // BANDGAP_H
