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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "LocalABO3EigenstressAux.h"

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", LocalABO3EigenstressAux);

InputParameters LocalABO3EigenstressAux::validParams() {
  InputParameters params = NodalPatchRecovery::validParams();
  params.addClassDescription("Local eigenstress calculation ");
  params.addRequiredCoupledVar("polar_x",
                               "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y",
                               "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i", "index_i >= 0 & index_i <= 2",
      "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j", "index_j >= 0 & index_j <= 2",
      "The index j of ij for the tensor to output (0, 1, 2)");
  return params;
}

LocalABO3EigenstressAux::LocalABO3EigenstressAux(
    const InputParameters &parameters)
    : NodalPatchRecovery(parameters), _i(getParam<unsigned int>("index_i")),
      _j(getParam<unsigned int>("index_j")), _polar_x(coupledValue("polar_x")),
      _polar_y(coupledValue("polar_y")), _polar_z(coupledValue("polar_z")),
      _C11(getMaterialProperty<Real>("C11")),
      _C12(getMaterialProperty<Real>("C12")),
      _C44(getMaterialProperty<Real>("C44")),
      _Q11(getMaterialProperty<Real>("Q11")),
      _Q12(getMaterialProperty<Real>("Q12")),
      _Q44(getMaterialProperty<Real>("Q44")) {}

Real LocalABO3EigenstressAux::computeValue() {

  RankTwoTensor eigenstress_tensor;

  eigenstress_tensor(0, 0) =
      _C11[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q11[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q11[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q11[_qp] +
      2.0 * _C12[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q12[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q12[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q12[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q12[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q12[_qp];

  eigenstress_tensor(1, 1) =
      _C12[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q11[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q11[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q11[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q12[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q12[_qp] +
      2.0 * _C12[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q12[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q12[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q12[_qp];

  eigenstress_tensor(2, 2) =
      _C12[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q11[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q11[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q11[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q12[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_x[_qp]) * _Q12[_qp] +
      _C11[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q12[_qp] +
      _C12[_qp] * Utility::pow<2>(_polar_y[_qp]) * _Q12[_qp] +
      2.0 * _C12[_qp] * Utility::pow<2>(_polar_z[_qp]) * _Q12[_qp];

  eigenstress_tensor(0, 1) = eigenstress_tensor(1, 0) =
      4.0 * _C44[_qp] * _polar_x[_qp] * _polar_y[_qp] * _Q44[_qp];

  eigenstress_tensor(1, 2) = eigenstress_tensor(2, 1) =
      4.0 * _C44[_qp] * _polar_y[_qp] * _polar_z[_qp] * _Q44[_qp];

  eigenstress_tensor(0, 2) = eigenstress_tensor(2, 0) =
      4.0 * _C44[_qp] * _polar_x[_qp] * _polar_z[_qp] * _Q44[_qp];

  if (_i == 0 & _j == 0)
    return eigenstress_tensor(0, 0);
  else if (_i == 1 & _j == 1)
    return eigenstress_tensor(1, 1);
  else if (_i == 2 & _j == 2)
    return eigenstress_tensor(2, 2);
  else if (_i == 0 & _j == 1)
    return eigenstress_tensor(0, 1);
  else if (_i == 1 & _j == 2)
    return eigenstress_tensor(1, 2);
  else if (_i == 0 & _j == 2)
    return eigenstress_tensor(0, 2);
  else
    return 0.0;
}
