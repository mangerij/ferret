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

#include "MasterZLSpinTransferTorqueCartLLG.h"

registerMooseObject("FerretApp", MasterZLSpinTransferTorqueCartLLG);

InputParameters
MasterZLSpinTransferTorqueCartLLG::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Zhang-Li spin-transfer torque for Cartesian LLG magnetization.");
  params.addRequiredParam<unsigned int>("component", "Magnetization component: 0=x, 1=y, 2=z.");
  params.addParam<unsigned int>("current_direction", 0, "Spatial direction of current: 0=x, 1=y, 2=z.");
  params.addRequiredCoupledVar("mag_x", "x component of magnetization.");
  params.addRequiredCoupledVar("mag_y", "y component of magnetization.");
  params.addRequiredCoupledVar("mag_z", "z component of magnetization.");
  params.addParam<Real>("u", 0.0, "Signed spin-drift velocity in mesh-length units per time unit.");
  params.addParam<Real>("xi", 0.0, "Non-adiabaticity parameter.");
  return params;
}

MasterZLSpinTransferTorqueCartLLG::MasterZLSpinTransferTorqueCartLLG(
    const InputParameters & parameters)
  : Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _current_direction(getParam<unsigned int>("current_direction")),
    _mag_x_var(coupled("mag_x")),
    _mag_y_var(coupled("mag_y")),
    _mag_z_var(coupled("mag_z")),
    _mag_x(coupledValue("mag_x")),
    _mag_y(coupledValue("mag_y")),
    _mag_z(coupledValue("mag_z")),
    _mag_x_grad(coupledGradient("mag_x")),
    _mag_y_grad(coupledGradient("mag_y")),
    _mag_z_grad(coupledGradient("mag_z")),
    _alpha(getMaterialProperty<Real>("alpha")),
    _u(getParam<Real>("u")),
    _xi(getParam<Real>("xi"))
{
  if (_component > 2)
    mooseError("component must be 0, 1, or 2.");

  if (_current_direction > 2)
    mooseError("current_direction must be 0, 1, or 2.");
}

Real
MasterZLSpinTransferTorqueCartLLG::computeQpResidual()
{
  return -torqueComponent(_component) * _test[_i][_qp];
}

Real
MasterZLSpinTransferTorqueCartLLG::computeQpJacobian()
{
  return -torqueDerivative(_component, _component) * _test[_i][_qp];
}

Real
MasterZLSpinTransferTorqueCartLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _mag_x_var)
    return -torqueDerivative(_component, 0) * _test[_i][_qp];

  if (jvar == _mag_y_var)
    return -torqueDerivative(_component, 1) * _test[_i][_qp];

  if (jvar == _mag_z_var)
    return -torqueDerivative(_component, 2) * _test[_i][_qp];

  return 0.0;
}

Real
MasterZLSpinTransferTorqueCartLLG::torqueComponent(unsigned int component) const
{
  const Real m[3] = {_mag_x[_qp], _mag_y[_qp], _mag_z[_qp]};
  const Real p[3] = {_mag_x_grad[_qp](_current_direction),
                     _mag_y_grad[_qp](_current_direction),
                     _mag_z_grad[_qp](_current_direction)};

  const Real m2 = m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
  const Real mdp = m[0] * p[0] + m[1] * p[1] + m[2] * p[2];

  const Real adv[3] = {p[0] * m2 - m[0] * mdp,
                       p[1] * m2 - m[1] * mdp,
                       p[2] * m2 - m[2] * mdp};

  const Real cross[3] = {m[1] * p[2] - m[2] * p[1],
                         m[2] * p[0] - m[0] * p[2],
                         m[0] * p[1] - m[1] * p[0]};

  const Real alpha = _alpha[_qp];
  const Real denominator = 1.0 + alpha * alpha;

  return (-_u * (1.0 + alpha * _xi) * adv[component] +
          _u * (_xi - alpha) * cross[component]) /
         denominator;
}

Real
MasterZLSpinTransferTorqueCartLLG::torqueDerivative(unsigned int component,
                                                  unsigned int variable_component) const
{
  const Real m[3] = {_mag_x[_qp], _mag_y[_qp], _mag_z[_qp]};
  const Real p[3] = {_mag_x_grad[_qp](_current_direction),
                     _mag_y_grad[_qp](_current_direction),
                     _mag_z_grad[_qp](_current_direction)};

  const Real m2 = m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
  const Real mdp = m[0] * p[0] + m[1] * p[1] + m[2] * p[2];
  const Real phi = _phi[_j][_qp];
  const Real grad_phi = _grad_phi[_j][_qp](_current_direction);

  const Real delta = component == variable_component ? 1.0 : 0.0;

  const Real d_adv_dm = 2.0 * p[component] * m[variable_component] -
                        delta * mdp - m[component] * p[variable_component];
  const Real d_adv_dp = delta * m2 - m[component] * m[variable_component];

  Real d_cross_dm = 0.0;
  Real d_cross_dp = 0.0;

  for (unsigned int k = 0; k < 3; ++k)
  {
    d_cross_dm += leviCivita(component, variable_component, k) * p[k];
    d_cross_dp += leviCivita(component, k, variable_component) * m[k];
  }

  const Real d_adv = d_adv_dm * phi + d_adv_dp * grad_phi;
  const Real d_cross = d_cross_dm * phi + d_cross_dp * grad_phi;

  const Real alpha = _alpha[_qp];
  const Real denominator = 1.0 + alpha * alpha;

  return (-_u * (1.0 + alpha * _xi) * d_adv + _u * (_xi - alpha) * d_cross) /
         denominator;
}

int
MasterZLSpinTransferTorqueCartLLG::leviCivita(unsigned int i, unsigned int j, unsigned int k) const
{
  if (i == j || j == k || i == k)
    return 0;

  if ((i == 0 && j == 1 && k == 2) || (i == 1 && j == 2 && k == 0) ||
      (i == 2 && j == 0 && k == 1))
    return 1;

  return -1;
}
