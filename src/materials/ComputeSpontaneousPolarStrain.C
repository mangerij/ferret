/*
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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ComputeSpontaneousPolarStrain.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeSpontaneousPolarStrain);

InputParameters ComputeSpontaneousPolarStrain::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Compute the spontaneous polar contribution to the strain.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

ComputeSpontaneousPolarStrain::ComputeSpontaneousPolarStrain(const InputParameters & parameters) :
    ComputeEigenstrainBase(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _Q11(getMaterialProperty<Real>("Q11")),
  _Q12(getMaterialProperty<Real>("Q12")),
  _Q44(getMaterialProperty<Real>("Q44")),
  _vals(6),
 _polar_strain()
{
}

void
ComputeSpontaneousPolarStrain::computeQpEigenstrain()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  _vals[0] = _Q11[_qp]*w(0)*w(0) + _Q12[_qp]*(w(1)*w(1)+w(2)*w(2)); // eps_1
  _vals[1] = _Q11[_qp]*w(1)*w(1) + _Q12[_qp]*(w(0)*w(0)+w(2)*w(2)); //eps_2
  _vals[2] = _Q11[_qp]*w(2)*w(2) + _Q12[_qp]*(w(0)*w(0)+w(1)*w(1)); //eps_3
  _vals[3] = _Q44[_qp]*(w(1)*w(2)); //eps_4   23
  _vals[4] = _Q44[_qp]*(w(0)*w(2)); //eps_5   13
  _vals[5] = _Q44[_qp]*(w(0)*w(1)); //eps_6   12
  _polar_strain.fillFromInputVector(_vals);
  _eigenstrain[_qp] = _polar_strain;
}
