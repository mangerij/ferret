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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ElasticEnergy.h"
#include "ComputeEigenstrain.h"

registerMooseObject("FerretApp", ElasticEnergy);

template<>
InputParameters validParams<ElasticEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the elastic energy density. Note this file also exists in tensor mechanics.");
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
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<" Stress-divergence equation:                                              "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"       ∇·[σ]  = 0                                                         "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";
}

Real
ElasticEnergy::computeQpIntegral()
{
  Real scaling = _len_scale*_len_scale*_len_scale*_strain_scale;
  return scaling*0.5*_stress[_qp].doubleContraction(_elastic_strain[_qp]);
}
