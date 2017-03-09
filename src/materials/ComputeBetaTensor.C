/****************************************************************/
/* Computes a rank 2 adjustment to indicatrix                   */
/****************************************************************/

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

