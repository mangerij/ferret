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

#ifndef ROTOBULKENERGYEIGHTH_H
#define ROTOBULKENERGYEIGHTH_H

#include "ElementIntegralPostprocessor.h"

class RotoBulkEnergyEighth;

template<>
InputParameters validParams<RotoBulkEnergyEighth>();

class RotoBulkEnergyEighth : public ElementIntegralPostprocessor
{
public:
  RotoBulkEnergyEighth(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _antiferrodis_A_x;
  const VariableValue& _antiferrodis_A_y;
  const VariableValue& _antiferrodis_A_z;
  const MaterialProperty<Real> & _beta1;
  const MaterialProperty<Real> & _beta11;
  const MaterialProperty<Real> & _beta12;
  const MaterialProperty<Real> & _beta111;
  const MaterialProperty<Real> & _beta112;
  const MaterialProperty<Real> & _beta123;
  const MaterialProperty<Real> & _beta1111;
  const MaterialProperty<Real> & _beta1112;
  const MaterialProperty<Real> & _beta1122;
  const MaterialProperty<Real> & _beta1123;

};

#endif
