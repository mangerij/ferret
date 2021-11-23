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

#ifndef BFOEIGENSTRAIN_H
#define BFOEIGENSTRAIN_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class BFOEigenstrain : public AuxKernel
{
public:
  BFOEigenstrain(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~BFOEigenstrain() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_j;
  const unsigned int _index_k;

  const VariableValue &_antiferrodis_A_x;
  const VariableValue &_antiferrodis_A_y;
  const VariableValue &_antiferrodis_A_z;
  const VariableValue &_polar_x;
  const VariableValue &_polar_y;
  const VariableValue &_polar_z;
  const Real _C11;
  const Real _C12;
  const Real _C44;
  const Real _Q11;
  const Real _Q12;
  const Real _Q44;
  const Real _R11;
  const Real _R12;
  const Real _R44;

  const MaterialProperty<RankTwoTensor> & _strain; //elastic?

};

#endif // BFOEIGENSTRAIN_H
