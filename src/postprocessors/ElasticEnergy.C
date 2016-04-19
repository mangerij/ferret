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
  params.addParam<Real>("strain_scale", 1.0, "the strain_scale");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

ElasticEnergy::ElasticEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  _strain_scale(getParam<Real>("strain_scale")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
ElasticEnergy::computeQpIntegral()
{
  return 0.5 * std::pow(_len_scale, 3.0) * _strain_scale * (
      _total_strain[_qp](0, 0) * _stress[_qp](0, 0)
    + _total_strain[_qp](0, 1) * _stress[_qp](0, 0)
    + _total_strain[_qp](0, 2) * _stress[_qp](0, 0)
    + _total_strain[_qp](1, 0) * _stress[_qp](0, 0)
    + _total_strain[_qp](1, 1) * _stress[_qp](0, 0)
    + _total_strain[_qp](1, 2) * _stress[_qp](0, 0)
    + _total_strain[_qp](2, 0) * _stress[_qp](0, 0)
    + _total_strain[_qp](2, 1) * _stress[_qp](0, 0)
    + _total_strain[_qp](2, 2) * _stress[_qp](0, 0)
    + _total_strain[_qp](0, 0) * _stress[_qp](0, 1)
    + _total_strain[_qp](0, 1) * _stress[_qp](0, 1)
    + _total_strain[_qp](0, 2) * _stress[_qp](0, 1)
    + _total_strain[_qp](1, 0) * _stress[_qp](0, 1)
    + _total_strain[_qp](1, 1) * _stress[_qp](0, 1)
    + _total_strain[_qp](1, 2) * _stress[_qp](0, 1)
    + _total_strain[_qp](2, 0) * _stress[_qp](0, 1)
    + _total_strain[_qp](2, 1) * _stress[_qp](0, 1)
    + _total_strain[_qp](2, 2) * _stress[_qp](0, 1)
    + _total_strain[_qp](0, 0) * _stress[_qp](0, 2)
    + _total_strain[_qp](0, 1) * _stress[_qp](0, 2)
    + _total_strain[_qp](0, 2) * _stress[_qp](0, 2)
    + _total_strain[_qp](1, 0) * _stress[_qp](0, 2)
    + _total_strain[_qp](1, 1) * _stress[_qp](0, 2)
    + _total_strain[_qp](1, 2) * _stress[_qp](0, 2)
    + _total_strain[_qp](2, 0) * _stress[_qp](0, 2)
    + _total_strain[_qp](2, 1) * _stress[_qp](0, 2)
    + _total_strain[_qp](2, 2) * _stress[_qp](0, 2)
    + _total_strain[_qp](0, 0) * _stress[_qp](1, 0)
    + _total_strain[_qp](0, 1) * _stress[_qp](1, 0)
    + _total_strain[_qp](0, 2) * _stress[_qp](1, 0)
    + _total_strain[_qp](1, 0) * _stress[_qp](1, 0)
    + _total_strain[_qp](1, 1) * _stress[_qp](1, 0)
    + _total_strain[_qp](1, 2) * _stress[_qp](1, 0)
    + _total_strain[_qp](2, 0) * _stress[_qp](1, 0)
    + _total_strain[_qp](2, 1) * _stress[_qp](1, 0)
    + _total_strain[_qp](2, 2) * _stress[_qp](1, 0)
    + _total_strain[_qp](0, 0) * _stress[_qp](1, 1)
    + _total_strain[_qp](0, 1) * _stress[_qp](1, 1)
    + _total_strain[_qp](0, 2) * _stress[_qp](1, 1)
    + _total_strain[_qp](1, 0) * _stress[_qp](1, 1)
    + _total_strain[_qp](1, 1) * _stress[_qp](1, 1)
    + _total_strain[_qp](1, 2) * _stress[_qp](1, 1)
    + _total_strain[_qp](2, 0) * _stress[_qp](1, 1)
    + _total_strain[_qp](2, 1) * _stress[_qp](1, 1)
    + _total_strain[_qp](2, 2) * _stress[_qp](1, 1)
    + _total_strain[_qp](0, 0) * _stress[_qp](1, 2)
    + _total_strain[_qp](0, 1) * _stress[_qp](1, 2)
    + _total_strain[_qp](0, 2) * _stress[_qp](1, 2)
    + _total_strain[_qp](1, 0) * _stress[_qp](1, 2)
    + _total_strain[_qp](1, 1) * _stress[_qp](1, 2)
    + _total_strain[_qp](1, 2) * _stress[_qp](1, 2)
    + _total_strain[_qp](2, 0) * _stress[_qp](1, 2)
    + _total_strain[_qp](2, 1) * _stress[_qp](1, 2)
    + _total_strain[_qp](2, 2) * _stress[_qp](1, 2)
    + _total_strain[_qp](0, 0) * _stress[_qp](2, 0)
    + _total_strain[_qp](0, 1) * _stress[_qp](2, 0)
    + _total_strain[_qp](0, 2) * _stress[_qp](2, 0)
    + _total_strain[_qp](1, 0) * _stress[_qp](2, 0)
    + _total_strain[_qp](1, 1) * _stress[_qp](2, 0)
    + _total_strain[_qp](1, 2) * _stress[_qp](2, 0)
    + _total_strain[_qp](2, 0) * _stress[_qp](2, 0)
    + _total_strain[_qp](2, 1) * _stress[_qp](2, 0)
    + _total_strain[_qp](2, 2) * _stress[_qp](2, 0)
    + _total_strain[_qp](0, 0) * _stress[_qp](2, 1)
    + _total_strain[_qp](0, 1) * _stress[_qp](2, 1)
    + _total_strain[_qp](0, 2) * _stress[_qp](2, 1)
    + _total_strain[_qp](1, 0) * _stress[_qp](2, 1)
    + _total_strain[_qp](1, 1) * _stress[_qp](2, 1)
    + _total_strain[_qp](1, 2) * _stress[_qp](2, 1)
    + _total_strain[_qp](2, 0) * _stress[_qp](2, 1)
    + _total_strain[_qp](2, 1) * _stress[_qp](2, 1)
    + _total_strain[_qp](2, 2) * _stress[_qp](2, 1)
    + _total_strain[_qp](0, 0) * _stress[_qp](2, 2)
    + _total_strain[_qp](0, 1) * _stress[_qp](2, 2)
    + _total_strain[_qp](0, 2) * _stress[_qp](2, 2)
    + _total_strain[_qp](1, 0) * _stress[_qp](2, 2)
    + _total_strain[_qp](1, 1) * _stress[_qp](2, 2)
    + _total_strain[_qp](1, 2) * _stress[_qp](2, 2)
    + _total_strain[_qp](2, 0) * _stress[_qp](2, 2)
    + _total_strain[_qp](2, 1) * _stress[_qp](2, 2)
    + _total_strain[_qp](2, 2) * _stress[_qp](2, 2)
  );
}
