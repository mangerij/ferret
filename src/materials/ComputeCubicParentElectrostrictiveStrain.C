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

#include "ComputeCubicParentElectrostrictiveStrain.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeCubicParentElectrostrictiveStrain);

InputParameters ComputeCubicParentElectrostrictiveStrain::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Compute the spontaneous polar contribution to the strain. It also declares its derivatives w.r.t. the coupled polarization components, so StressDivergenceTensors assembles the u<-P Jacobian block when given eigenstrain_names = <this eigenstrain_name> and coupled_variables = <the polarization components>.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

ComputeCubicParentElectrostrictiveStrain::ComputeCubicParentElectrostrictiveStrain(const InputParameters & parameters) :
    DerivativeMaterialInterface<ComputeEigenstrainBase>(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _Q11(getMaterialProperty<Real>("Q11")),
  _Q12(getMaterialProperty<Real>("Q12")),
  _Q44(getMaterialProperty<Real>("Q44")),
  _vals(6),
  _cubic_parent_electrostrictive_strain(),
  _deigenstrain(3, nullptr)
{
  const std::vector<std::string> comps = {"polar_x", "polar_y", "polar_z"};
  for (unsigned int k = 0; k < 3; ++k)
    if (isCoupled(comps[k]))
      _deigenstrain[k] = &declarePropertyDerivative<RankTwoTensor>(_eigenstrain_name, coupledName(comps[k]));
}

void
ComputeCubicParentElectrostrictiveStrain::computeQpEigenstrain()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  _vals[0] = _Q11[_qp]*w(0)*w(0) + _Q12[_qp]*(w(1)*w(1)+w(2)*w(2)); // eps_1
  _vals[1] = _Q11[_qp]*w(1)*w(1) + _Q12[_qp]*(w(0)*w(0)+w(2)*w(2)); //eps_2
  _vals[2] = _Q11[_qp]*w(2)*w(2) + _Q12[_qp]*(w(0)*w(0)+w(1)*w(1)); //eps_3
  _vals[3] = _Q44[_qp]*(w(1)*w(2)); //eps_4   23
  _vals[4] = _Q44[_qp]*(w(0)*w(2)); //eps_5   13
  _vals[5] = _Q44[_qp]*(w(0)*w(1)); //eps_6   12
  _cubic_parent_electrostrictive_strain.fillFromInputVector(_vals);
  _eigenstrain[_qp] = _cubic_parent_electrostrictive_strain;

  // d eps0_ij / d P_k: the diagonal from Q11/Q12, the (k,j) shear pair from Q44
  for (unsigned int k = 0; k < 3; ++k)
    if (_deigenstrain[k])
    {
      RankTwoTensor & D = (*_deigenstrain[k])[_qp];
      D.zero();
      for (unsigned int j = 0; j < 3; ++j)
        if (j == k)
          D(j, j) = 2.0 * _Q11[_qp] * w(k);
        else
        {
          D(j, j) = 2.0 * _Q12[_qp] * w(k);
          D(k, j) = D(j, k) = _Q44[_qp] * w(j);
        }
    }
}
