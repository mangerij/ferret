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

#include "MagneticExchangeEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticExchangeEnergy);

template<>
InputParameters validParams<MagneticExchangeEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the magnetic exchange energy density.");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization");
  return params;
}

MagneticExchangeEnergy::MagneticExchangeEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _mag_x_grad(coupledGradient("mag_x")),
  _mag_y_grad(coupledGradient("mag_y")),
  _mag_z_grad(coupledGradient("mag_z")),
  _Ae(getMaterialProperty<Real>("Ae")),
  _mu0(getMaterialProperty<Real>("mu0"))
{
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"Selecting:                                                                "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<" Landau-Liftshitz-Bloch equations for evolution of the magnetic system    "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<" dmk/dt = - γ' (mk × δF/δmk - α [mk × (mk × δF/δmk] + αL [m^4 - m^2] m_k) "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  //TODO: later can rework this in the following way: postprocessors will print energetic contributions and a "blank" kernel will print the LGD/LLG/coupled terms
  //      can also use for elastic and electrostatic coupling.
}

Real
MagneticExchangeEnergy::computeQpIntegral()
{
  return _mu0[_qp]*((_Ae[_qp])*(Utility::pow<2>(_mag_x_grad[_qp](0))+Utility::pow<2>(_mag_x_grad[_qp](1))+Utility::pow<2>(_mag_x_grad[_qp](2))+Utility::pow<2>(_mag_y_grad[_qp](0))+Utility::pow<2>(_mag_y_grad[_qp](1))+Utility::pow<2>(_mag_y_grad[_qp](2))+Utility::pow<2>(_mag_z_grad[_qp](0))+Utility::pow<2>(_mag_z_grad[_qp](1))+Utility::pow<2>(_mag_z_grad[_qp](2))));
}
