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

#include "AFDWallEnergy.h"

registerMooseObject("FerretApp", AFDWallEnergy);

InputParameters AFDWallEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the computational volume of the free energy density due to afd vector field gradients"
                             "corresponding to gradients in the AFD field.");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the afd vector field");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

AFDWallEnergy::AFDWallEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _antiphase_A_x_grad(coupledGradient("antiphase_A_x")),
  _antiphase_A_y_grad(coupledGradient("antiphase_A_y")),
  _antiphase_A_z_grad(coupledGradient("antiphase_A_z")),
  _H110(getMaterialProperty<Real>("H110")),
  _H11(getMaterialProperty<Real>("H11_H110")),
  _H12(getMaterialProperty<Real>("H12_H110")),
  _H44(getMaterialProperty<Real>("H44_H110")),
  _H44P(getMaterialProperty<Real>("H44P_H110")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
AFDWallEnergy::computeQpIntegral()
{
  return _energy_scale*(_H110[_qp]*(0.5*_H11[_qp]*(Utility::pow<2>(_antiphase_A_x_grad[_qp](0))+Utility::pow<2>(_antiphase_A_y_grad[_qp](1))+Utility::pow<2>(_antiphase_A_z_grad[_qp](2)))+
    _H12[_qp]*(_antiphase_A_x_grad[_qp](0)*_antiphase_A_y_grad[_qp](1)+_antiphase_A_y_grad[_qp](1)*_antiphase_A_z_grad[_qp](2)+_antiphase_A_x_grad[_qp](0)*_antiphase_A_z_grad[_qp](2))+
    0.5*_H44[_qp]*(pow(_antiphase_A_x_grad[_qp](1)+_antiphase_A_y_grad[_qp](0),2)+Utility::pow<2>(_antiphase_A_y_grad[_qp](2)+_antiphase_A_z_grad[_qp](1))+Utility::pow<2>(_antiphase_A_x_grad[_qp](2)+_antiphase_A_z_grad[_qp](0))+
	  0.5*_H44P[_qp]*(Utility::pow<2>(_antiphase_A_x_grad[_qp](1)-_antiphase_A_y_grad[_qp](0))+Utility::pow<2>(_antiphase_A_y_grad[_qp](2)-_antiphase_A_z_grad[_qp](1))+Utility::pow<2>(_antiphase_A_x_grad[_qp](2)-_antiphase_A_z_grad[_qp](0))))));
}
