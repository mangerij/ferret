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

#include "MagneticExchangeEnergyUS.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticExchangeEnergyUS);

InputParameters MagneticExchangeEnergyUS::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the magnetic exchange energy density.");
  params.addRequiredCoupledVar("azimuthal_ph", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_th", "The polar component of the constrained magnetic vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and aJ");
  return params;
}

MagneticExchangeEnergyUS::MagneticExchangeEnergyUS(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _azimuthal_ph(coupledValue("azimuthal_ph")),
  _polar_th(coupledValue("polar_th")),
  _polar_th_grad(coupledGradient("polar_th")),
  _azimuthal_ph_grad(coupledGradient("azimuthal_ph")),
  _Ae(getMaterialProperty<Real>("Ae")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _energy_scale(getParam<Real>("energy_scale"))
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
MagneticExchangeEnergyUS::computeQpIntegral()
{
  return _energy_scale*((_Ae[_qp]*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + 2.0*(Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2))) - (Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)))*std::cos(2.0*_polar_th[_qp])))/2.);
}
