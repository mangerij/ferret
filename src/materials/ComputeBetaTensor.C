/**
 * @file   ComputeBetaTensor.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate photoelastic change to the refractive index
 * \Delta \epsilon_{ij} = - n_i^2 n_j^2 \Delta B_{ij}
 *
 * where \Delta B_{ij} = p_{ijkl} \varepsilon_{kl}
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#include "ComputeBetaTensor.h"

template<>
InputParameters validParams<ComputeBetaTensor>()
{
  InputParameters params = validParams<ComputeBetaTensorBase>();
  params.addClassDescription("Compute the adjustments to the indicatrix (beta tensor).");
  return params;
}

ComputeBetaTensor::ComputeBetaTensor(const InputParameters & parameters) :
    ComputeBetaTensorBase(parameters),
    _strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
    _photostrictive_tensor(getMaterialProperty<RankFourTensor>("photostrictive_tensor"))
{
}

void
ComputeBetaTensor::computeQpBetaTensor()
{
  for (unsigned int a = 0; a < 3; ++a)
    for (unsigned int b = 0; b < 3; ++b)
    {
      Real sum = 0;
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
        {
          sum += _photostrictive_tensor[_qp](a, b, i,j) * _strain[_qp](i,j);
        }
    _beta_tensor[_qp](a, b) = sum;
    }
}

