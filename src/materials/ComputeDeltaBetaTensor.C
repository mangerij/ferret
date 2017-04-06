/**
 * @file   ComputeDeltaBetaTensor.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate photoelastic change to the refractive index
 *
 * where \Delta (1/n^2) = \Delta B_{ij} = p_{ijkl} \varepsilon_{kl}
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#include "ComputeDeltaBetaTensor.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<ComputeDeltaBetaTensor>()
{
  InputParameters params = validParams<ComputeDeltaBetaTensorBase>();
  params.addClassDescription("Compute the adjustments to the indicatrix (beta tensor).");
  return params;
}

ComputeDeltaBetaTensor::ComputeDeltaBetaTensor(const InputParameters & parameters) :
    ComputeDeltaBetaTensorBase(parameters),
    _strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
    _photostrictive_tensor(getMaterialProperty<RankFourTensor>("photostrictive_tensor"))
{
}

void
ComputeDeltaBetaTensor::computeQpDeltaBetaTensor()
{
  Real sum = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
        {
          sum += _photostrictive_tensor[_qp](i, j, k, l) * _strain[_qp](k,l);
        }
    _delta_beta_tensor[_qp](i, j) = sum;
    }
    //Moose::out << "\n b"; std::cout << a; Moose::out << " = "; std::cout << _delta_beta_tensor[_qp](0, a);
}


