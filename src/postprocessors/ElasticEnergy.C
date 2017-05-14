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

#include "ElasticEnergy.h"
#include "ComputeEigenstrain.h"

template<>
InputParameters validParams<ElasticEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addParam<Real>("strain_scale", 1.0, "the strain_scale");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

ElasticEnergy::ElasticEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _elastic_strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  _strain_scale(getParam<Real>("strain_scale")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
ElasticEnergy::computeQpIntegral()
{
  Real scaling = _len_scale*_len_scale*_len_scale*_strain_scale;

  return scaling*0.5*_stress[_qp].doubleContraction(_elastic_strain[_qp]);

}
