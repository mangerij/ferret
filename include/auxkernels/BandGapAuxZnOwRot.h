#ifndef BANDGAPAUXZNOWROT_H
#define BANDGAPAUXZNOWROT_H

#include "AuxKernel.h"
#include "Material.h"
#include "TensorMechanicsMaterial.h" //may not need this
#include "RankTwoTensor.h"

//Forward declarations
class BandGapAuxZnOwRot;

template<>
InputParameters validParams<BandGapAuxZnOwRot>();


class BandGapAuxZnOwRot : public AuxKernel
{
public:
  BandGapAuxZnOwRot(const InputParameters & parameters);

  virtual ~BandGapAuxZnOwRot() {}

protected:
  virtual Real computeValue();

private:
  const MaterialProperty<RankTwoTensor> & _strain;
  const Real _du, _db, _E0, _Rb, _nu;
  RealVectorValue _Euler_angles;
};

#endif // BANDGAP_H
