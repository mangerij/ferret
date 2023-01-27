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

#include "HarmonicFieldAux.h"

registerMooseObject("FerretApp", HarmonicFieldAux);

InputParameters HarmonicFieldAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates a harmonic field");
  params.addRequiredParam<Real>("amplitude", "amplitude of the field");
  params.addParam<Real>("correction", 1.0, "correction (for rotations)");
  params.addRequiredParam<Real>("frequency", "frequency of the field");
  params.addParam<Real>("tshift", 0.0, "shift of the field");
  params.addParam<Real>("ton", 0.0, "switch on time of the field");
  params.addParam<Real>("toff", 0.0, "switch off time of the field");
  return params;
}


HarmonicFieldAux::HarmonicFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _amplitude(getParam<Real>("amplitude")),
   _correction(getParam<Real>("correction")),
   _frequency(getParam<Real>("frequency")),
   _tshift(getParam<Real>("tshift")),
   _ton(getParam<Real>("ton")),
   _toff(getParam<Real>("toff"))
{
}

Real
HarmonicFieldAux::computeValue()
{
  /*Moose::out << "\n time = "; std::cout << _t;
  Moose::out << "\n field turns on at "; std::cout << _ton;
  Moose::out << "\n field turns off at "; std::cout << _toff;*/
  if (_t > _ton && _t < _toff)
  {
    return _amplitude*_correction*std::sin(_frequency*(_t+_tshift));
  }
  else
    return 0.0;
}


