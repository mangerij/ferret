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

#ifndef BULKENERGYDENSITY_H
#define BULKENERGYDENSITY_H

#include "AuxKernel.h"

class BulkEnergyDensity;

template<>
InputParameters validParams<BulkEnergyDensity>();

class BulkEnergyDensity : public AuxKernel
{
public:
  BulkEnergyDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const MaterialProperty<Real> & _alpha1;
  const MaterialProperty<Real> & _alpha11;
  const MaterialProperty<Real> & _alpha12;
  const MaterialProperty<Real> & _alpha111;
  const MaterialProperty<Real> & _alpha112;
  const MaterialProperty<Real> & _alpha123;
  const MaterialProperty<Real> & _alpha1111;
  const MaterialProperty<Real> & _alpha1112;
  const MaterialProperty<Real> & _alpha1122;
  const MaterialProperty<Real> & _alpha1123;
};

#endif // BULKENERGYDENSITY_H
