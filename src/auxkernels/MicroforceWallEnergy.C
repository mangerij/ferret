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

#include "MicroforceWallEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MicroforceWallEnergy);

template<>
InputParameters validParams<MicroforceWallEnergy>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

MicroforceWallEnergy::MicroforceWallEnergy(const InputParameters & parameters) :
  AuxKernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_i_grad((_component==0)? coupledGradient("polar_x") :(_component==1)? coupledGradient("polar_y"): coupledGradient("polar_z")),
  _polar_j_grad((_component==0)? coupledGradient("polar_y"): (_component==1)? coupledGradient("polar_z"): coupledGradient("polar_x")),
  _polar_k_grad((_component==0)? coupledGradient("polar_z"): (_component==1)? coupledGradient("polar_x"): coupledGradient("polar_y")),
  _ii(_component),
  _jj((_component==0)? 1 : (_component==1)? 2: 0),
  _kk((_component==0)? 2 : (_component==1)? 0: 1),
  _G110(getMaterialProperty<Real>("G110")),
  _G11(getMaterialProperty<Real>("G11_G110")),
  _G12(getMaterialProperty<Real>("G12_G110")),
  _G44(getMaterialProperty<Real>("G44_G110")),
  _G44P(getMaterialProperty<Real>("G44P_G110"))
{}

Real
MicroforceWallEnergy::computeValue()
{
//important: note that this return is consistent with the Kernel and NOT the AuxKernel calculating the density/postprocessor (factor of 2
// also, this is merely an approximation since I set _grad_test[_i][_qp](_kk) -> 1.0 for now. TODO: calc actual func derivative without integration by parts
  Real Rwall = 0.0;

  Rwall += (_G11[_qp] * _polar_i_grad[_qp](_ii) * 1.0 +
    _G12[_qp] * (_polar_j_grad[_qp](_jj) + _polar_k_grad[_qp](_kk)) * 1.0 +
    _G44[_qp] * (_polar_i_grad[_qp](_jj) + _polar_j_grad[_qp](_ii)) * 1.0 + _G44[_qp] * (_polar_i_grad[_qp](_kk)+_polar_k_grad[_qp](_ii)) * 1.0 +
	   _G44P[_qp] * (_polar_i_grad[_qp](_jj) - _polar_j_grad[_qp](_ii)) * 1.0 + _G44P[_qp] * (_polar_i_grad[_qp](_kk) - _polar_k_grad[_qp](_ii)) * 1.0);
  ///  Moose::out << "\n R_wall-"; std::cout << _component << " = " << Rwall;
  return _G110[_qp]*Rwall;
}
