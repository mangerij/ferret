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
  for (unsigned int a = 0; a < 3; ++a)
    {
      Real sum = 0;
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
        {
          sum += _photostrictive_tensor[_qp](a, a, i, j) * _strain[_qp](i,j); //This needs to match what we have in Mathematica (Nye notation).
          //Moose::out << "\n b"; std::cout << a << i << j; Moose::out << " = "; std::cout << _photostrictive_tensor[_qp](a, a, i, j) * _strain[_qp](i,j);
        }
      _delta_beta_tensor[_qp](0, a) = sum;
      //Moose::out << "\n b"; std::cout << a; Moose::out << " = "; std::cout << _delta_beta_tensor[_qp](0, a);
    }
    Real sum1 = 0;
    Real sum2 = 0;
    Real sum3 = 0;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
      {
        sum1 += _photostrictive_tensor[_qp](1, 2, i, j) * _strain[_qp](i,j);
        sum2 += _photostrictive_tensor[_qp](2, 0, i, j) * _strain[_qp](i,j);
        sum3 += _photostrictive_tensor[_qp](0, 1, i, j) * _strain[_qp](i,j);
      }
  _delta_beta_tensor[_qp](0,3) = sum1;
  _delta_beta_tensor[_qp](0,4) = sum2;
  _delta_beta_tensor[_qp](0,5) = sum3;
}


