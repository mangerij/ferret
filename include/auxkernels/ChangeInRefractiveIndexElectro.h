/**
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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef CHANGEINREFRACTIVEINDEXELECTRO_H
#define CHANGEINREFRACTIVEINDEXELECTRO_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class ChangeInRefractiveIndexElectro;

template<>
InputParameters validParams<ChangeInRefractiveIndexElectro>();


class ChangeInRefractiveIndexElectro : public AuxKernel
{
public:
  ChangeInRefractiveIndexElectro(const InputParameters & parameters);

  virtual ~ChangeInRefractiveIndexElectro() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_i;
  const unsigned int _index_j;
  const unsigned int _index_k;
  const unsigned int _index_l;
  const MaterialProperty<RankTwoTensor> & _indicatrix;
  const MaterialProperty<RankTwoTensor> & _delta_indicatrix_electro;
};

#endif // CHANGEINREFRACTIVEINDEXELECTRO_H
