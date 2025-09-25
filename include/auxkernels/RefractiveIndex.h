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

#ifndef REFRACTIVEINDEX_H
#define REFRACTIVEINDEX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class RefractiveIndex : public AuxKernel
{
public:
  RefractiveIndex(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~RefractiveIndex() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_j;
  const unsigned int _index_k;
  bool _electro;
  bool _elasto;
  bool _polar;
  const MaterialProperty<RankTwoTensor> & _indicatrix;
  const VariableValue & _var1;
  const VariableValue & _var2;
  const VariableValue & _var3;
};

#endif // REFRACTIVEINDEX_H
