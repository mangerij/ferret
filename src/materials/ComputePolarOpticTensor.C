/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

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
  _elastooptic_tensor(getMaterialProperty<RankFourTensor>("elastooptic_tensor")),
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
          sum += _elastooptic_tensor[_qp](i, j, k, l) * _electrostrictive_coefficients[_qp](k, l, m, n) * w(m) * w(n);
          }
        }
    _delta_PO_tensor[_qp](i, j) = sum;
    }
}


