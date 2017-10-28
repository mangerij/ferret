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

#include "ComputePolarOpticGCoeffTensor.h"
#include "RankTwoTensor.h"

template<>
InputParameters validParams<ComputePolarOpticGCoeffTensor>()
{
  InputParameters params = validParams<ComputePolarOpticGCoeffTensorBase>();
  params.addClassDescription("Compute the adjustments to the indicatrix due to the polar-optic effect with gijkl coefficients.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

ComputePolarOpticGCoeffTensor::ComputePolarOpticGCoeffTensor(const InputParameters & parameters) :
   ComputePolarOpticGCoeffTensorBase(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _gcoefficient_tensor(getMaterialProperty<RankFourTensor>("gcoefficient_tensor"))
{
}

void
ComputePolarOpticGCoeffTensor::computeQpPolarOpticGCoeffTensor()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  Real sum = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      for (unsigned int m = 0; m < 3; ++m)
        for (unsigned int n = 0; n < 3; ++n)
        {
          sum += _gcoefficient_tensor[_qp](i, j, m, n) * w(m) * w(n);
        }
    _delta_gPO_tensor[_qp](i, j) = sum;
    }
}


