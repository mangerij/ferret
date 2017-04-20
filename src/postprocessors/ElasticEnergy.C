/**
 * @file   ElasticEnergy.C
 * @author J. Mangeri <mangerij@anl.gov>
 * @modified A. Jookisaari <andrea.jokisaari@northwestern.edu>
 *
 * @brief This is a energy postprocessor that tracks the elastic energy
 *        which is equivalent to one half the contraction of stress on 
 *        strain or 1/2 * \sigma_{ij} \varepsilon_{ij}.
 */

#include "ElasticEnergy.h"
#include "ComputeEigenstrain.h"

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
  _elastic_strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  _strain_scale(getParam<Real>("strain_scale")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
ElasticEnergy::computeQpIntegral()
{
  Real scaling = _len_scale*_len_scale*_len_scale*_strain_scale;

  return scaling*0.5*_stress[_qp].doubleContraction(_elastic_strain[_qp]);

}
