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

#include "MasterInteractionCartLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MasterInteractionCartLLG);

InputParameters MasterInteractionCartLLG::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the Rij contribution (due to energy -M*H), assuming H = - div*Phi.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", 0.0, "The external magnetic potential variable");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addParam<Real>("g0", 1.0, "electron gyromagnetic factor");
  params.addParam<Real>("Hscale", 1.0, "scaling factor for effective fields");
  return params;
}

MasterInteractionCartLLG::MasterInteractionCartLLG(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_H_int_var(coupled("potential_H_int")),
   _potential_H_ext_var(coupled("potential_H_ext")),
   _potential_H_int(coupledValue("potential_H_int")),
   _potential_H_ext(coupledValue("potential_H_ext")),
   _potential_H_int_grad(coupledGradient("potential_H_int")),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _alpha(getMaterialProperty<Real>("alpha")),
   _g0(getParam<Real>("g0")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _Hscale(getParam<Real>("Hscale"))
{
}


Real
MasterInteractionCartLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return -((_g0/_Hscale)*(_potential_H_int_grad[_qp](2)*_mag_y[_qp] - _potential_H_int_grad[_qp](1)*_mag_z[_qp] + _alpha[_qp]*1.0*(_potential_H_int_grad[_qp](1)*_mag_x[_qp]*_mag_y[_qp] + _potential_H_int_grad[_qp](2)*_mag_x[_qp]*_mag_z[_qp] - _potential_H_int_grad[_qp](0)*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 1)
  {
    return -((_g0/_Hscale)*(-(_potential_H_int_grad[_qp](2)*_mag_x[_qp]) + _potential_H_int_grad[_qp](0)*_mag_z[_qp] + _alpha[_qp]*1.0*(_mag_y[_qp]*(_potential_H_int_grad[_qp](0)*_mag_x[_qp] + _potential_H_int_grad[_qp](2)*_mag_z[_qp]) - _potential_H_int_grad[_qp](1)*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 2)
  {
    return -((_g0/_Hscale)*(_potential_H_int_grad[_qp](1)*_mag_x[_qp] - _potential_H_int_grad[_qp](0)*_mag_y[_qp] + _alpha[_qp]*1.0*(-(_potential_H_int_grad[_qp](2)*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))) + (_potential_H_int_grad[_qp](0)*_mag_x[_qp] + _potential_H_int_grad[_qp](1)*_mag_y[_qp])*_mag_z[_qp]))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else
    return 0.0;
}

Real
MasterInteractionCartLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return -(_alpha[_qp]*(_g0/_Hscale)*1.0*(_potential_H_int_grad[_qp](1)*_mag_y[_qp] + _potential_H_int_grad[_qp](2)*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 1)
  {
    return -(_alpha[_qp]*(_g0/_Hscale)*1.0*(_potential_H_int_grad[_qp](0)*_mag_x[_qp] + _potential_H_int_grad[_qp](2)*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 2)
  {
    return -(_alpha[_qp]*(_g0/_Hscale)*1.0*(_potential_H_int_grad[_qp](0)*_mag_x[_qp] + _potential_H_int_grad[_qp](1)*_mag_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else
    return 0.0;
}

Real
MasterInteractionCartLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return -((_g0/_Hscale)*(_potential_H_int_grad[_qp](2) + _alpha[_qp]*1.0*(_potential_H_int_grad[_qp](1)*_mag_x[_qp] - 2*_potential_H_int_grad[_qp](0)*_mag_y[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return -((_g0/_Hscale)*(-_potential_H_int_grad[_qp](1) + _alpha[_qp]*1.0*(_potential_H_int_grad[_qp](2)*_mag_x[_qp] - 2*_potential_H_int_grad[_qp](0)*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_int_var)
    {
      return -((_g0/_Hscale)*(_grad_phi[_j][_qp](2)*_mag_y[_qp] - _grad_phi[_j][_qp](1)*_mag_z[_qp] + _alpha[_qp]*1.0*(_grad_phi[_j][_qp](1)*_mag_x[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](2)*_mag_x[_qp]*_mag_z[_qp] - _grad_phi[_j][_qp](0)*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
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
    if (jvar == _mag_x_var)
    {
      return -((_g0/_Hscale)*(-_potential_H_int_grad[_qp](2) + _alpha[_qp]*1.0*(-2*_potential_H_int_grad[_qp](1)*_mag_x[_qp] + _potential_H_int_grad[_qp](0)*_mag_y[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return -((_g0/_Hscale)*(_potential_H_int_grad[_qp](0) + _alpha[_qp]*1.0*(_potential_H_int_grad[_qp](2)*_mag_y[_qp] - 2*_potential_H_int_grad[_qp](1)*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_int_var)
    {
      return -((_g0/_Hscale)*(-(_grad_phi[_j][_qp](2)*_mag_x[_qp]) + _grad_phi[_j][_qp](0)*_mag_z[_qp] + _alpha[_qp]*1.0*(_mag_y[_qp]*(_grad_phi[_j][_qp](0)*_mag_x[_qp] + _grad_phi[_j][_qp](2)*_mag_z[_qp]) - _grad_phi[_j][_qp](1)*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
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
  else if (_component == 2)
  {
    if (jvar == _mag_x_var)
    {
      return -((_g0/_Hscale)*(_potential_H_int_grad[_qp](1) + _alpha[_qp]*1.0*(-2*_potential_H_int_grad[_qp](2)*_mag_x[_qp] + _potential_H_int_grad[_qp](0)*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
      return -((_g0/_Hscale)*(-_potential_H_int_grad[_qp](0) + _alpha[_qp]*1.0*(-2*_potential_H_int_grad[_qp](2)*_mag_y[_qp] + _potential_H_int_grad[_qp](1)*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _potential_H_int_var)
    {
      return -((_g0/_Hscale)*(_grad_phi[_j][_qp](1)*_mag_x[_qp] - _grad_phi[_j][_qp](0)*_mag_y[_qp] + _alpha[_qp]*1.0*(-(_grad_phi[_j][_qp](2)*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))) + (_grad_phi[_j][_qp](0)*_mag_x[_qp] + _grad_phi[_j][_qp](1)*_mag_y[_qp])*_mag_z[_qp]))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
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
