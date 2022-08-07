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

#include "MasterInteractionUSLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MasterInteractionUSLLG);

InputParameters MasterInteractionUSLLG::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution - M$*$H in the total energy, assuming H = - div * potential.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for polar theta, 1 for azimuthal phi");
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", 0.0, "The external magnetic potential variable");
  params.addRequiredCoupledVar("polar_th", "The polar angle of the constrained magnetic vector");
  params.addRequiredCoupledVar("azimuthal_ph", "The azimuthal angle of the constrained magnetic vector");
  return params;
}

MasterInteractionUSLLG::MasterInteractionUSLLG(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_H_int_var(coupled("potential_H_int")),
   _potential_H_ext_var(coupled("potential_H_ext")),
   _potential_H_int(coupledValue("potential_H_int")),
   _potential_H_ext(coupledValue("potential_H_ext")),
   _potential_H_int_grad(coupledGradient("potential_H_int")),
   _polar_th_var(coupled("polar_th")),
   _azimuthal_ph_var(coupled("azimuthal_ph")),
   _polar_th(coupledValue("polar_th")),
   _azimuthal_ph(coupledValue("azimuthal_ph")),
   _alpha(getMaterialProperty<Real>("alpha")),
   _g0(getMaterialProperty<Real>("g0")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}


Real
MasterInteractionUSLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return -((_g0[_qp]*_test[_i][_qp]*(std::cos(_azimuthal_ph[_qp])*(_potential_H_int_grad[_qp](1) + _alpha[_qp]*_potential_H_int_grad[_qp](0)*std::cos(_polar_th[_qp])) + (-_potential_H_int_grad[_qp](0) + _alpha[_qp]*_potential_H_int_grad[_qp](1)*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp]) - _alpha[_qp]*_potential_H_int_grad[_qp](2)*std::sin(_polar_th[_qp])))/(1 + Utility::pow<2>(_alpha[_qp])));
  }
  else if (_component == 1)
  {
    return -((_g0[_qp]*_test[_i][_qp]*(_potential_H_int_grad[_qp](2) + std::cos(_azimuthal_ph[_qp])*(_alpha[_qp]*_potential_H_int_grad[_qp](1) + _potential_H_int_grad[_qp](0)*std::cos(_polar_th[_qp]))*(1.0/std::sin(_polar_th[_qp])) - (_alpha[_qp]*_potential_H_int_grad[_qp](0) + _potential_H_int_grad[_qp](1)*std::cos(_polar_th[_qp]))*(1.0/std::sin(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp])))/
     (1 + Utility::pow<2>(_alpha[_qp])));
  }
  else
    return 0.0;
}

Real
MasterInteractionUSLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_alpha[_qp]*_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(_potential_H_int_grad[_qp](2)*std::cos(_polar_th[_qp]) + (_potential_H_int_grad[_qp](0)*std::cos(_azimuthal_ph[_qp]) + _potential_H_int_grad[_qp](1)*std::sin(_azimuthal_ph[_qp]))*std::sin(_polar_th[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 1)
  {
    return (_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(1.0/std::sin(_polar_th[_qp]))*(std::cos(_azimuthal_ph[_qp])*(_alpha[_qp]*_potential_H_int_grad[_qp](0) + _potential_H_int_grad[_qp](1)*std::cos(_polar_th[_qp])) + (_alpha[_qp]*_potential_H_int_grad[_qp](1) + _potential_H_int_grad[_qp](0)*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else
    return 0.0;
}

Real
MasterInteractionUSLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuthal_ph_var)
    {
      return (_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(std::cos(_azimuthal_ph[_qp])*(_potential_H_int_grad[_qp](0) - _alpha[_qp]*_potential_H_int_grad[_qp](1)*std::cos(_polar_th[_qp])) + (_potential_H_int_grad[_qp](1) + _alpha[_qp]*_potential_H_int_grad[_qp](0)*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_int_var)
    {
      return (_g0[_qp]*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))*(std::cos(_azimuthal_ph[_qp])*(_potential_H_int_grad[_qp](0) + _alpha[_qp]*_potential_H_int_grad[_qp](1)*std::cos(_polar_th[_qp])) - (_potential_H_int_grad[_qp](1) + _alpha[_qp]*_potential_H_int_grad[_qp](0)*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp])))/(1 + Utility::pow<2>(_alpha[_qp]));
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
    if (jvar == _polar_th_var)
    {
      return -((_g0[_qp]*_test[_i][_qp]*(_grad_phi[_j][_qp](2) + std::cos(_azimuthal_ph[_qp])*(_alpha[_qp]*_grad_phi[_j][_qp](1) + _grad_phi[_j][_qp](0)*std::cos(_polar_th[_qp]))*(1.0/std::sin(_polar_th[_qp])) - (_alpha[_qp]*_grad_phi[_j][_qp](0) + _grad_phi[_j][_qp](1)*std::cos(_polar_th[_qp]))*(1.0/std::sin(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp])))/
     (1 + Utility::pow<2>(_alpha[_qp])));
    }
    else if (jvar == _potential_H_int_var)
    {
      return -((_g0[_qp]*_test[_i][_qp]*(std::cos(_azimuthal_ph[_qp])*(_grad_phi[_j][_qp](1) + _alpha[_qp]*_grad_phi[_j][_qp](0)*std::cos(_polar_th[_qp])) + (-_grad_phi[_j][_qp](0) + _alpha[_qp]*_grad_phi[_j][_qp](1)*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](2)*std::sin(_polar_th[_qp])))/(1 + Utility::pow<2>(_alpha[_qp])));
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
