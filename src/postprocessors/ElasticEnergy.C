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
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

ElasticEnergy::ElasticEnergy(const std::string & name, InputParameters parameters) :
  ElementIntegralPostprocessor(name, parameters),
  _elastic_strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
ElasticEnergy::computeQpIntegral()
{
  return -0.5 * std::pow(_len_scale, 3.0) * (
      _elastic_strain[_qp](0, 0) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](0, 0)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](0, 1)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](0, 2)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](1, 0)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](1, 1)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](1, 2)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](2, 0)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](2, 1)
    + _elastic_strain[_qp](0, 0) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](0, 1) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](0, 2) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](1, 0) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](1, 1) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](1, 2) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](2, 0) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](2, 1) * _stress[_qp](2, 2)
    + _elastic_strain[_qp](2, 2) * _stress[_qp](2, 2)
  );
}
