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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ComputeSpontaneousRotostrictiveStrain.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeSpontaneousRotostrictiveStrain);

InputParameters ComputeSpontaneousRotostrictiveStrain::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Compute the rotostrictive contributuon to the spontaneous strain. It also declares its derivatives w.r.t. the coupled antiphase tilt components, so StressDivergenceTensors assembles the u<-A Jacobian block when given eigenstrain_names = <this eigenstrain_name> and coupled_variables = <the antiphase tilt components>.");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt order parameter");
  params.addCoupledVar("antiphase_A_y", 0.0, "The y component of the antiphase tilt order parameter");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt order parameter");
  return params;
}

ComputeSpontaneousRotostrictiveStrain::ComputeSpontaneousRotostrictiveStrain(const InputParameters & parameters) :
    DerivativeMaterialInterface<ComputeEigenstrainBase>(parameters),
  _antiphase_A_x(coupledValue("antiphase_A_x")),
  _antiphase_A_y(coupledValue("antiphase_A_y")),
  _antiphase_A_z(coupledValue("antiphase_A_z")),
  _R11(getMaterialProperty<Real>("R11")),
  _R12(getMaterialProperty<Real>("R12")),
  _R44(getMaterialProperty<Real>("R44")),
  _vals(6),
 _rotostrictive_strain(),
  _deigenstrain(3, nullptr)
{
  const std::vector<std::string> comps = {"antiphase_A_x", "antiphase_A_y", "antiphase_A_z"};
  for (unsigned int k = 0; k < 3; ++k)
    if (isCoupled(comps[k]))
      _deigenstrain[k] = &declarePropertyDerivative<RankTwoTensor>(_eigenstrain_name, coupledName(comps[k]));
}

void
ComputeSpontaneousRotostrictiveStrain::computeQpEigenstrain()
{
  RealVectorValue w(_antiphase_A_x[_qp], _antiphase_A_y[_qp], _antiphase_A_z[_qp]);

  _vals[0] = _R11[_qp]*w(0)*w(0) + _R12[_qp]*(w(1)*w(1)+w(2)*w(2)); // eps_1
  _vals[1] = _R11[_qp]*w(1)*w(1) + _R12[_qp]*(w(0)*w(0)+w(2)*w(2)); //eps_2
  _vals[2] = _R11[_qp]*w(2)*w(2) + _R12[_qp]*(w(0)*w(0)+w(1)*w(1)); //eps_3
  _vals[3] = _R44[_qp]*(w(1)*w(2)); //eps_4   23
  _vals[4] = _R44[_qp]*(w(0)*w(2)); //eps_5   13
  _vals[5] = _R44[_qp]*(w(0)*w(1)); //eps_6   12
  _rotostrictive_strain.fillFromInputVector(_vals);
  _eigenstrain[_qp] = _rotostrictive_strain;

  // d eps0_ij / d A_k: the diagonal from R11/R12, the (k,j) shear pair from R44
  for (unsigned int k = 0; k < 3; ++k)
    if (_deigenstrain[k])
    {
      RankTwoTensor & D = (*_deigenstrain[k])[_qp];
      D.zero();
      for (unsigned int j = 0; j < 3; ++j)
        if (j == k)
          D(j, j) = 2.0 * _R11[_qp] * w(k);
        else
        {
          D(j, j) = 2.0 * _R12[_qp] * w(k);
          D(k, j) = D(j, k) = _R44[_qp] * w(j);
        }
    }
}
