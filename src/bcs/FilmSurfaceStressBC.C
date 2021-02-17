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

#include "FilmSurfaceStressBC.h"

registerMooseObject("FerretApp", FilmSurfaceStressBC);

template<>
InputParameters validParams<FilmSurfaceStressBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addClassDescription(
      "stress free surface condition, testing, only for z direction");
  params.addRequiredParam<int>("component","Which component(0 for x, 1 for y, 2 for z) in traction is used");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

FilmSurfaceStressBC::FilmSurfaceStressBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _component(getParam<int>("component")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _polar_x(coupledValue("polar_x")),
    _polar_y(coupledValue("polar_y")),
    _polar_z(coupledValue("polar_z")),
    _C11(getMaterialProperty<Real>("C11")),
    _C12(getMaterialProperty<Real>("C12")),
    _C44(getMaterialProperty<Real>("C44")),
    _Q11(getMaterialProperty<Real>("Q11")),
    _Q12(getMaterialProperty<Real>("Q12")),
    _Q44(getMaterialProperty<Real>("Q44"))
{
}

Real
FilmSurfaceStressBC::computeQpResidual()
{
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

  return _test[_i][_qp] * (_stress[_qp](_component,2)+eigenstress_tensor(_component,2));
}

Real
FilmSurfaceStressBC::computeQpJacobian()
{
  return 0;
}
