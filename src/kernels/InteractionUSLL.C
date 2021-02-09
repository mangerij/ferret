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

#include "InteractionUSLL.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", InteractionUSLL);

template<>
InputParameters validParams<InteractionUSLL>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution - M$*$H in the total energy, assuming H = - div * potential.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for polar, 1 for azimuth)");
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", 0.0, "The external magnetic potential variable");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  return params;
}

InteractionUSLL::InteractionUSLL(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_H_int_var(coupled("potential_H_int")),
   _potential_H_ext_var(coupled("potential_H_ext")),
   _potential_H_int(coupledValue("potential_H_int")),
   _potential_H_ext(coupledValue("potential_H_ext")),
   _potential_H_int_grad(coupledGradient("potential_H_int")),
   _azimuth_phi_var(coupled("azimuth_phi")),
   _polar_theta_var(coupled("polar_theta")),
   _azimuth_phi(coupledValue("azimuth_phi")),
   _polar_theta(coupledValue("polar_theta")),
   _alpha(getMaterialProperty<Real>("alpha")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _g0(getMaterialProperty<Real>("g0"))
{
}


Real
InteractionUSLL::computeQpResidual()
{
  if (_component == 0)
  {
    return (_g0[_qp]*_test[_i][_qp]*std::sin(_polar_theta[_qp])*(-(std::cos(_azimuth_phi[_qp])*(_potential_H_int_grad[_qp](1) + _alpha[_qp]*_potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp]))) + (_potential_H_int_grad[_qp](0) - _alpha[_qp]*_potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _alpha[_qp]*_potential_H_int_grad[_qp](2)*std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 1)
  {
    return (_g0[_qp]*_test[_i][_qp]*(std::cos(_azimuth_phi[_qp])*(_alpha[_qp]*_potential_H_int_grad[_qp](1) - _potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp])) - (_alpha[_qp]*_potential_H_int_grad[_qp](0) + _potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _potential_H_int_grad[_qp](2)*std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else
    return 0.0;
}

Real
InteractionUSLL::computeQpJacobian()
{
  if (_component == 0)
  {
    return -((_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(std::cos(_polar_theta[_qp])*(_potential_H_int_grad[_qp](1)*std::cos(_azimuth_phi[_qp]) - _potential_H_int_grad[_qp](0)*std::sin(_azimuth_phi[_qp])) + _alpha[_qp]*std::cos(2*_polar_theta[_qp])*(_potential_H_int_grad[_qp](0)*std::cos(_azimuth_phi[_qp]) + _potential_H_int_grad[_qp](1)*std::sin(_azimuth_phi[_qp])) - _alpha[_qp]*_potential_H_int_grad[_qp](2)*std::sin(2*_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha[_qp])));
  }
  else if (_component == 1)
  {
    return -((_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(std::cos(_azimuth_phi[_qp])*(_alpha[_qp]*_potential_H_int_grad[_qp](0) + _potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp])) + (_alpha[_qp]*_potential_H_int_grad[_qp](1) - _potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])))/(1 + Utility::pow<2>(_alpha[_qp])));
  }
  else
    return 0.0;
}

Real
InteractionUSLL::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return (_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(std::cos(_azimuth_phi[_qp])*(_potential_H_int_grad[_qp](0) - _alpha[_qp]*_potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp])) + (_potential_H_int_grad[_qp](1) + _alpha[_qp]*_potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]))/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_int_var)
    {
      return (_g0[_qp]*_test[_i][_qp]*std::sin(_polar_theta[_qp])*(-(std::cos(_azimuth_phi[_qp])*(_grad_phi[_j][_qp](1) + _alpha[_qp]*_grad_phi[_j][_qp](0)*std::cos(_polar_theta[_qp]))) + (_grad_phi[_j][_qp](0) - _alpha[_qp]*_grad_phi[_j][_qp](1)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _alpha[_qp]*_grad_phi[_j][_qp](2)*std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_ext_var)
    {
      return 0.0;
    }
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (jvar == _polar_theta_var)
    {
      return (_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(_potential_H_int_grad[_qp](2)*std::cos(_polar_theta[_qp]) + (_potential_H_int_grad[_qp](0)*std::cos(_azimuth_phi[_qp]) + _potential_H_int_grad[_qp](1)*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_int_var)
    {
      return (_g0[_qp]*_test[_i][_qp]*(std::cos(_azimuth_phi[_qp])*(_alpha[_qp]*_grad_phi[_j][_qp](1) - _grad_phi[_j][_qp](0)*std::cos(_polar_theta[_qp])) - (_alpha[_qp]*_grad_phi[_j][_qp](0) + _grad_phi[_j][_qp](1)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _grad_phi[_j][_qp](2)*std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_ext_var)
    {
      return 0.0;
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
