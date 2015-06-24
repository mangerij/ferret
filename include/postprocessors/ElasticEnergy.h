/**
 * @file   WallEnergy.h
 * @author J. Mangeri <mangerij@anl.gov>
 *
 */

#ifndef ELASTICENERGY_H
#define ELASTICENERGY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class ElasticEnergy;

template<>
InputParameters validParams<ElasticEnergy>();


class ElasticEnergy : public ElementIntegralPostprocessor
{
public:
  ElasticEnergy(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpIntegral();

private:
  const MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _stress;

};

#endif
