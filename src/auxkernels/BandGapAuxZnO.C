/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "BandGapAuxZnO.h"

template<>

InputParameters validParams<BandGapAuxZnO>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("relaxed_energy", 0.0,"relaxed energy");
  params.addParam<Real>("biaxial_strain_rate", 0.0, "uniaxial strain rate");
  params.addParam<Real>("uniaxial_strain_rate", 0.0, "biaxial strain rate");
  params.addParam<Real>("biaxial_relaxation_coeff", 0.0, "biaxial relaxation coeff");
  params.addParam<Real>("poisson_ratio", 0.0, "Poisson ratio");
  return params;
}


BandGapAuxZnO::BandGapAuxZnO(const InputParameters & parameters) :
  AuxKernel(parameters),
   _strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
   _du(getParam<Real>("biaxial_strain_rate")),
   _db(getParam<Real>("uniaxial_strain_rate")),
   _E0(getParam<Real>("relaxed_energy")),
   _Rb(getParam<Real>("biaxial_relaxation_coeff")),
   _nu(getParam<Real>("poisson_ratio"))
{
}

Real
BandGapAuxZnO::computeValue()

{
    return _E0 + (1/(1-_Rb))*((_db+_du*_Rb)*0.5*(_strain[_qp](1,1)+_strain[_qp](0,0))+(_du+_nu*_db)*_strain[_qp](2,2));
}


