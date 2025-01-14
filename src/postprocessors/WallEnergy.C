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

#include "WallEnergy.h"

registerMooseObject("FerretApp", WallEnergy);

InputParameters WallEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the Ginzburg term.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

WallEnergy::WallEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G110(getMaterialProperty<Real>("G110")),
  _G11(getMaterialProperty<Real>("G11_G110")),
  _G12(getMaterialProperty<Real>("G12_G110")),
  _G44(getMaterialProperty<Real>("G44_G110")),
  _G44P(getMaterialProperty<Real>("G44P_G110")),
  _energy_scale(getParam<Real>("energy_scale"))
{
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<" time-dependent coupled Landau-Ginzburg equations for evolution of the    "<<"\n";
  std::cout<<" ferroelectric system:                                                    "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"       dPk/dt = -Γ δF/δPk                                                 "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";
  //TODO: later can rework this in the following way: postprocessors will print energetic contributions and a "blank" postprocessor will print the LGD/LLG/coupled terms
  //      can also use for elastic and electrostatic coupling.
}

Real
WallEnergy::computeQpIntegral()
{
  return _energy_scale*_G110[_qp]*(0.5*_G11[_qp]*
                   (
                    pow(_polar_x_grad[_qp](0),2)+pow(_polar_y_grad[_qp](1),2)+pow(_polar_z_grad[_qp](2),2)
                   )
+              _G12[_qp]*(
                    _polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)+_polar_y_grad[_qp](1)*_polar_z_grad[_qp](2)+_polar_x_grad[_qp](0)*_polar_z_grad[_qp](2)
                   )
+          0.5*_G44[_qp]*(
                    pow(_polar_x_grad[_qp](1)+_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)+_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)+_polar_z_grad[_qp](0),2)
                   )
+
	  0.5*_G44P[_qp]*(pow(_polar_x_grad[_qp](1)-_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)-_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)-_polar_z_grad[_qp](0),2)));
}
