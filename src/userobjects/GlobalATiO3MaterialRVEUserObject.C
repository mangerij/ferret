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

#include "GlobalATiO3MaterialRVEUserObject.h"

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", GlobalATiO3MaterialRVEUserObject);

InputParameters
GlobalATiO3MaterialRVEUserObject::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription(
      "Global Strain UserObject to provide Residual and diagonal Jacobian entry");
  params.addParam<std::vector<Real>>("applied_stress_tensor",
                                     "Vector of values defining the constant applied stress "
                                     "to add, in order 11, 22, 33, 23, 13, 12");
  params.addParam<std::string>("base_name", "Material properties base name");
  params.addCoupledVar("displacements", "The name of the displacement variables");
  params.set<ExecFlagEnum>("execute_on") = EXEC_LINEAR;
  params.addCoupledVar("polar_x", 0.0, "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

GlobalATiO3MaterialRVEUserObject::GlobalATiO3MaterialRVEUserObject(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _dstress_dstrain(getMaterialProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")),
    _dim(_mesh.dimension()),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _periodic_dir(),
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
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  for (unsigned int dir = 0; dir < _dim; ++dir)
  {
    _periodic_dir(dir) = _mesh.isTranslatedPeriodic(_disp_var[0], dir);

    for (unsigned int i = 1; i < _ndisp; ++i)
      if (_mesh.isTranslatedPeriodic(_disp_var[i], dir) != _periodic_dir(dir))
        mooseError("All the displacement components in a particular direction should have same "
                   "periodicity.");
  }

  if (isParamValid("applied_stress_tensor"))
    _applied_stress_tensor.fillFromInputVector(
        getParam<std::vector<Real>>("applied_stress_tensor"));
  else
    _applied_stress_tensor.zero();
}

void
GlobalATiO3MaterialRVEUserObject::initialize()
{
  _residual.zero();
  _jacobian.zero();
}

void
GlobalATiO3MaterialRVEUserObject::execute()
{
  computeAdditionalStress();

  for (unsigned int _qp = 0; _qp < _qrule->n_points(); _qp++)
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

    // residual, integral of stress components
    _residual += _JxW[_qp] * _coord[_qp] * (_stress[_qp] - _applied_stress_tensor + eigenstress_tensor);

    // diagonal jacobian, integral of elasticity tensor components
    _jacobian += _JxW[_qp] * _coord[_qp] * _dstress_dstrain[_qp];
  }
}

void
GlobalATiO3MaterialRVEUserObject::threadJoin(const UserObject & uo)
{
  const GlobalATiO3MaterialRVEUserObject & pstuo = static_cast<const GlobalATiO3MaterialRVEUserObject &>(uo);
  _residual += pstuo._residual;
  _jacobian += pstuo._jacobian;
}

void
GlobalATiO3MaterialRVEUserObject::finalize()
{
  std::vector<Real> residual(9);
  std::vector<Real> jacobian(81);

  std::copy(&_residual(0, 0), &_residual(0, 0) + 9, residual.begin());
  std::copy(&_jacobian(0, 0, 0, 0), &_jacobian(0, 0, 0, 0) + 81, jacobian.begin());

  gatherSum(residual);
  gatherSum(jacobian);

  std::copy(residual.begin(), residual.end(), &_residual(0, 0));
  std::copy(jacobian.begin(), jacobian.end(), &_jacobian(0, 0, 0, 0));
}

const RankTwoTensor &
GlobalATiO3MaterialRVEUserObject::getResidual() const
{
  return _residual;
}

const RankFourTensor &
GlobalATiO3MaterialRVEUserObject::getJacobian() const
{
  return _jacobian;
}

const VectorValue<bool> &
GlobalATiO3MaterialRVEUserObject::getPeriodicDirections() const
{
  return _periodic_dir;
}
