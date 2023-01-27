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

#ifndef REWORKEDREFRACTIVEINDEX_H
#define REWORKEDREFRACTIVEINDEX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class ReworkedRefractiveIndex : public AuxKernel
{
public:
  ReworkedRefractiveIndex(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~ReworkedRefractiveIndex() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _component;
  const MaterialProperty<Real> & _n1;
  const MaterialProperty<Real> & _n2;
  const MaterialProperty<Real> & _n3;
  const MaterialProperty<Real> & _n4;
  const MaterialProperty<Real> & _n5;
  const MaterialProperty<Real> & _n6;
  bool _electro;
  bool _elasto;
  bool _polar;
  const VariableValue & _var1;
  const VariableValue & _var2;
  const VariableValue & _var3;
};

#endif // REWORKEDREFRACTIVEINDEX_H
