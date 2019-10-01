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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "InteractionUSLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", InteractionUSLLG);

template<>
InputParameters validParams<InteractionUSLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution - M$*$H in the total energy, assuming H = - div * potential.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for polar, 1 for azimuth)");
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", 0.0, "The external magnetic potential variable");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("g0", "g0");
  return params;
}

InteractionUSLLG::InteractionUSLLG(const InputParameters & parameters)
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
   _alpha(getParam<Real>("alpha")),
   _g0(getParam<Real>("g0"))
{
}


Real
InteractionUSLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return -((_g0*_test[_i][_qp]*(std::cos(_azimuth_phi[_qp])*(_potential_H_int_grad[_qp](1) + _alpha*_potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp])) + (-_potential_H_int_grad[_qp](0) + _alpha*_potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) - _alpha*_potential_H_int_grad[_qp](2)*std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha)));
  }
  else if (_component == 1)
  {
    return (_g0*_test[_i][_qp]*(-_potential_H_int_grad[_qp](2) + std::cos(_azimuth_phi[_qp])*(-(_alpha*_potential_H_int_grad[_qp](1)) + _potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp]))*(1.0/std::sin(_polar_theta[_qp])) + _potential_H_int_grad[_qp](1)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _alpha*_potential_H_int_grad[_qp](0)*(1.0/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])))/(1 + Utility::pow<2>(_alpha));
  }
  else
    return 0.0;
}

Real
InteractionUSLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_alpha*_g0*_phi[_j][_qp]*_test[_i][_qp]*(_potential_H_int_grad[_qp](0)*std::cos(_azimuth_phi[_qp]) + _potential_H_int_grad[_qp](2)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + _potential_H_int_grad[_qp](1)*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]))/(1 + Utility::pow<2>(_alpha));
  }
  else if (_component == 1)
  {
    return (_g0*_phi[_j][_qp]*_test[_i][_qp]*(1.0/std::sin(_polar_theta[_qp]))*(std::cos(_azimuth_phi[_qp])*(_alpha*_potential_H_int_grad[_qp](0) + _potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp])) + (_alpha*_potential_H_int_grad[_qp](1) - _potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])))/(1 + Utility::pow<2>(_alpha));
  }
  else
    return 0.0;
}

Real
InteractionUSLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return (_g0*_phi[_j][_qp]*_test[_i][_qp]*(std::cos(_azimuth_phi[_qp])*(_potential_H_int_grad[_qp](0) - _alpha*_potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp])) + (_potential_H_int_grad[_qp](1) + _alpha*_potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])))/(1 + Utility::pow<2>(_alpha));
    }
    else if (jvar == _potential_H_int_var)
    {
      return -((_g0*_test[_i][_qp]*(std::cos(_azimuth_phi[_qp])*(_grad_phi[_j][_qp](1) + _alpha*_grad_phi[_j][_qp](0)*std::cos(_polar_theta[_qp])) + (-_grad_phi[_j][_qp](0) + _alpha*_grad_phi[_j][_qp](1)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) - _alpha*_grad_phi[_j][_qp](2)*std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha)));
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
      return -((_g0*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp])))*(std::cos(_azimuth_phi[_qp])*(_potential_H_int_grad[_qp](0) - _alpha*_potential_H_int_grad[_qp](1)*std::cos(_polar_theta[_qp])) + (_potential_H_int_grad[_qp](1) + _alpha*_potential_H_int_grad[_qp](0)*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])))/(1 + Utility::pow<2>(_alpha)));
    }
    else if (jvar == _potential_H_int_var)
    {
      return (_g0*_test[_i][_qp]*(-_grad_phi[_j][_qp](2) + std::cos(_azimuth_phi[_qp])*(-(_alpha*_grad_phi[_j][_qp](1)) + _grad_phi[_j][_qp](0)*std::cos(_polar_theta[_qp]))*(1.0/std::sin(_polar_theta[_qp])) + _grad_phi[_j][_qp](1)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _alpha*_grad_phi[_j][_qp](0)*(1.0/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])))/(1 + Utility::pow<2>(_alpha));
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
