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

#ifndef POLAROPTICCHANGEINREFRACTIVEINDEX_H
#define POLAROPTICCHANGEINREFRACTIVEINDEX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class PolarOpticChangeInRefractiveIndex : public AuxKernel
{
public:
  PolarOpticChangeInRefractiveIndex(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~PolarOpticChangeInRefractiveIndex() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _component;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _n1;
  const MaterialProperty<Real> & _n2;
  const MaterialProperty<Real> & _n3;
  const MaterialProperty<Real> & _n4;
  const MaterialProperty<Real> & _n5;
  const MaterialProperty<Real> & _n6;
  const MaterialProperty<Real> & _g1111;
  const MaterialProperty<Real> & _g1122;
  const MaterialProperty<Real> & _g1212;
};

#endif // POLAROPTICCHANGEINREFRACTIVEINDEX_H
