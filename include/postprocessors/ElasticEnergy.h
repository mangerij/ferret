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

#ifndef ELASTICENERGY_H
#define ELASTICENERGY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ElementIntegralPostprocessor.h"
#include "ComputeEigenstrain.h"

//Forward Declarations
class ElasticEnergy;

template<>
InputParameters validParams<ElasticEnergy>();


class ElasticEnergy : public ElementIntegralPostprocessor
{
public:
  ElasticEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

private:
  const MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _stress;
  const Real _strain_scale;
  const Real _len_scale;
};

#endif
