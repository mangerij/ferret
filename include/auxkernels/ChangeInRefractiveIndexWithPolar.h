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

#ifndef CHANGEINREFRACTIVEINDEXWITHPOLAR_H
#define CHANGEINREFRACTIVEINDEXWITHPOLAR_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class ChangeInRefractiveIndexWithPolar : public AuxKernel
{
public:
  ChangeInRefractiveIndexWithPolar(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~ChangeInRefractiveIndexWithPolar() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_i;
  const unsigned int _index_j;
  const unsigned int _index_k;
  const unsigned int _index_l;
  const MaterialProperty<RankTwoTensor> & _indicatrix;
  const MaterialProperty<RankTwoTensor> & _delta_indicatrix;
  const MaterialProperty<RankTwoTensor> & _delta_PO_tensor;
};

#endif // CHANGEINREFRACTIVEINDEXWITHPOLAR_H
