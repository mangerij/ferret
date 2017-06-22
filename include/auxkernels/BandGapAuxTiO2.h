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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

*/

#ifndef BANDGAPAUXTIO2_H
#define BANDGAPAUXTIO2_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class BandGapAuxTiO2;

template<>
InputParameters validParams<BandGapAuxTiO2>();


class BandGapAuxTiO2 : public AuxKernel
{
public:
  BandGapAuxTiO2(const InputParameters & parameters);

  virtual ~BandGapAuxTiO2() {}

protected:
  virtual Real computeValue();

private:
  const MaterialProperty<RankTwoTensor> & _stress;
  const Real _ba, _bc, _E0;

};

#endif // BANDGAP_H
