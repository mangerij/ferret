//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GlobalATOMaterialRVEUserObject.h"

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", GlobalATOMaterialRVEUserObject);

template <>
InputParameters
validParams<GlobalATOMaterialRVEUserObject>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription(
      "Global Strain UserObject to provide Residual and diagonal Jacobian entry");
  params.addParam<std::vector<Real>>("applied_stress_tensor",
                                     "Vector of values defining the constant applied stress "
                                     "to add, in order 11, 22, 33, 23, 13, 12");
  params.addParam<std::string>("base_name", "Material properties base name");
  params.addCoupledVar("displacements", "The name of the displacement variables");
  params.set<ExecFlagEnum>("execute_on") = EXEC_LINEAR;
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("C11", "the 11 component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C12", "the 12 component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C44", "the 44 component of elastic stiffness tensor");
  params.addRequiredParam<Real>("Q11", "the 11 component of electrostrictive coupling tensor");
  params.addRequiredParam<Real>("Q12", "the 12 component of electrostrictive coupling tensor");
  params.addRequiredParam<Real>("Q44", "the 44 component of electrostrictive coupling tensor");
  return params;
}

GlobalATOMaterialRVEUserObject::GlobalATOMaterialRVEUserObject(const InputParameters & parameters)
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
    _C11(getParam<Real>("C11")),
    _C12(getParam<Real>("C12")),
    _C44(getParam<Real>("C44")),
    _Q11(getParam<Real>("Q11")),
    _Q12(getParam<Real>("Q12")),
    _Q44(getParam<Real>("Q44"))
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
GlobalATOMaterialRVEUserObject::initialize()
{
  _residual.zero();
  _jacobian.zero();
}

void
GlobalATOMaterialRVEUserObject::execute()
{
  computeAdditionalStress();

  for (unsigned int _qp = 0; _qp < _qrule->n_points(); _qp++)
  {
   RankTwoTensor eigenstress_tensor;

   eigenstress_tensor(0, 0) =
        _C11 * Utility::pow<2>(_polar_x[_qp]) * _Q11 +
        _C12 * Utility::pow<2>(_polar_y[_qp]) * _Q11 +
        _C12 * Utility::pow<2>(_polar_z[_qp]) * _Q11 +
        2.0 * _C12 * Utility::pow<2>(_polar_x[_qp]) * _Q12 +
        _C11 * Utility::pow<2>(_polar_y[_qp]) * _Q12 +
        _C12 * Utility::pow<2>(_polar_y[_qp]) * _Q12 +
        _C11 * Utility::pow<2>(_polar_z[_qp]) * _Q12 +
        _C12 * Utility::pow<2>(_polar_z[_qp]) * _Q12;

    eigenstress_tensor(1, 1) =
        _C12 * Utility::pow<2>(_polar_x[_qp]) * _Q11 +
        _C11 * Utility::pow<2>(_polar_y[_qp]) * _Q11 +
        _C12 * Utility::pow<2>(_polar_z[_qp]) * _Q11 +
        _C11 * Utility::pow<2>(_polar_x[_qp]) * _Q12 +
        _C12 * Utility::pow<2>(_polar_x[_qp]) * _Q12 +
        2.0 * _C12 * Utility::pow<2>(_polar_y[_qp]) * _Q12 +
        _C11 * Utility::pow<2>(_polar_z[_qp]) * _Q12 +
        _C12 * Utility::pow<2>(_polar_z[_qp]) * _Q12;

    eigenstress_tensor(2, 2) =
        _C12 * Utility::pow<2>(_polar_x[_qp]) * _Q11 +
        _C12 * Utility::pow<2>(_polar_y[_qp]) * _Q11 +
        _C11 * Utility::pow<2>(_polar_z[_qp]) * _Q11 +
        _C11 * Utility::pow<2>(_polar_x[_qp]) * _Q12 +
        _C12 * Utility::pow<2>(_polar_x[_qp]) * _Q12 +
        _C11 * Utility::pow<2>(_polar_y[_qp]) * _Q12 +
        _C12 * Utility::pow<2>(_polar_y[_qp]) * _Q12 +
        2.0 * _C12 * Utility::pow<2>(_polar_z[_qp]) * _Q12;

    eigenstress_tensor(0, 1) = eigenstress_tensor(1, 0) =
        4.0 * _C44 * _polar_x[_qp] * _polar_y[_qp] * _Q44;

    eigenstress_tensor(1, 2) = eigenstress_tensor(2, 1) =
        4.0 * _C44 * _polar_y[_qp] * _polar_z[_qp] * _Q44;

    eigenstress_tensor(0, 2) = eigenstress_tensor(2, 0) =
        4.0 * _C44 * _polar_x[_qp] * _polar_z[_qp] * _Q44;

    // residual, integral of stress components
    _residual += _JxW[_qp] * _coord[_qp] * (_stress[_qp] - _applied_stress_tensor + eigenstress_tensor);

    // diagonal jacobian, integral of elasticity tensor components
    _jacobian += _JxW[_qp] * _coord[_qp] * _dstress_dstrain[_qp];
  }
}

void
GlobalATOMaterialRVEUserObject::threadJoin(const UserObject & uo)
{
  const GlobalATOMaterialRVEUserObject & pstuo = static_cast<const GlobalATOMaterialRVEUserObject &>(uo);
  _residual += pstuo._residual;
  _jacobian += pstuo._jacobian;
}

void
GlobalATOMaterialRVEUserObject::finalize()
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
GlobalATOMaterialRVEUserObject::getResidual() const
{
  return _residual;
}

const RankFourTensor &
GlobalATOMaterialRVEUserObject::getJacobian() const
{
  return _jacobian;
}

const VectorValue<bool> &
GlobalATOMaterialRVEUserObject::getPeriodicDirections() const
{
  return _periodic_dir;
}
