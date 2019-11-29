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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "Intensity.h"

registerMooseObject("FerretApp", Intensity);

template<>

InputParameters validParams<Intensity>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("ReE_x", "The (real) x component of the scattered electric field");
  params.addRequiredCoupledVar("ReE_y", "The (real) y component of the scattered electric field");
  params.addRequiredCoupledVar("ReE_z", "The (real) z component of the scattered electric field");
  params.addRequiredCoupledVar("ImagE_x", "The (imaginary) x component of the scattered electric field");
  params.addRequiredCoupledVar("ImagE_y", "The (imaginary) y component of the scattered electric field");
  params.addRequiredCoupledVar("ImagE_z", "The (imaginary) z component of the scattered electric field");
 // params.addRequiredCoupledVar("ReH_x", "The (real) x component of the scattered magnetic field");
 // params.addRequiredCoupledVar("ReH_y", "The (real) y component of the scattered magnetic field");
 // params.addRequiredCoupledVar("ReH_z", "The (real) z component of the scattered magnetic field");
 // params.addRequiredCoupledVar("ImagH_x", "The (imaginary) x component of the scattered magnetic field");
 // params.addRequiredCoupledVar("ImagH_y", "The (imaginary) y component of the scattered magnetic field");
 // params.addRequiredCoupledVar("ImagH_z", "The (imaginary) z component of the scattered magnetic field");
  return params;
}


Intensity::Intensity(const InputParameters & parameters) :
  AuxKernel(parameters),
   _ReE_x(coupledValue("ReE_x")),
   _ReE_y(coupledValue("ReE_y")),
   _ReE_z(coupledValue("ReE_z")),
   _ImagE_x(coupledValue("ImagE_x")),
   _ImagE_y(coupledValue("ImagE_y")),
   _ImagE_z(coupledValue("ImagE_z"))
  // _ReH_x(coupledValue("ReH_x")),
  // _ReH_y(coupledValue("ReH_y")),
  // _ReH_z(coupledValue("ReH_z")),
  // _ImagH_x(coupledValue("ImagH_x")),
  // _ImagH_y(coupledValue("ImagH_y")),
 //  _ImagH_z(coupledValue("ImagH_z"))
{
}

Real
Intensity::computeValue()
{
 return std::log(std::abs(_ReE_x[_qp] * _ImagE_x[_qp] + _ReE_y[_qp] * _ImagE_y[_qp] + _ReE_z[_qp] * _ImagE_z[_qp]));
 //   return (_ReE_x[_qp] * _ImagH_x[_qp] + _ReE_y[_qp] * _ImagH_y[_qp] + _ReE_z[_qp] * _ImagH_z[_qp]) + (_ImagE_x[_qp] * _ReH_x[_qp] + _ImagE_y[_qp] * _ReH_y[_qp] + _ImagE_z[_qp] * _ReH_z[_qp]);
}
