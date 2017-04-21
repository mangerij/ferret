/**
 * @file   ComputePolarOpticTensor.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate polar-optic change to the refractive index
 * \delta B_{ij} = p_{ijkl} Q_{klmn} P_m P_n
 *
 */


#include "ComputePolarOpticTensor.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<ComputePolarOpticTensor>()
{
  InputParameters params = validParams<ComputePolarOpticTensorBase>();
  params.addClassDescription("Compute the adjustments to the indicatrix (beta tensor) due to the polar-optic effect.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

ComputePolarOpticTensor::ComputePolarOpticTensor(const InputParameters & parameters) :
   ComputePolarOpticTensorBase(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _photostrictive_tensor(getMaterialProperty<RankFourTensor>("photostrictive_tensor")),
  _electrostrictive_coefficients(getMaterialProperty<RankFourTensor>("electrostrictive_coefficients"))
{
}

void
ComputePolarOpticTensor::computeQpPolarOpticTensor()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  Real sum = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
        {
        for (unsigned int m = 0; m < 3; ++m)
          for (unsigned int n = 0; n < 3; ++n)
          {
          sum += _photostrictive_tensor[_qp](i, j, k, l) * _electrostrictive_coefficients[_qp](k, l, m, n) * w(m) * w(n);
          }
        }
    _delta_PO_tensor[_qp](i, j) = sum;
    }
}


