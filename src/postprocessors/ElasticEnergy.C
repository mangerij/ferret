/**
 * @file   ElasticEnergy.C
 * @author J. Manger <mangerij@anl.gov>
 *
 */

#include "ElasticEnergy.h"

template<>
InputParameters validParams<ElasticEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
//  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

ElasticEnergy::ElasticEnergy(const std::string & name, InputParameters parameters) :
  ElementIntegralPostprocessor(name, parameters),
  _stress( getMaterialProperty<RankTwoTensor>("stress") ),
  _elastic_strain(getMaterialProperty<RankTwoTensor>("elastic_strain"))
{}

Real
ElasticEnergy::computeQpIntegral()
{
  return 0.5 * _stress[_qp].doubleContraction(_elastic_strain[_qp]);;
}
