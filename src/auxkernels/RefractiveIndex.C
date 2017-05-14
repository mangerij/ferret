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

#include "RefractiveIndex.h"

template<>

InputParameters validParams<RefractiveIndex>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("component", "component");
  params.addRequiredParam<Real>("n_a", "alpha refractive index");
  params.addRequiredParam<Real>("n_b", "beta refractive index");
  params.addRequiredParam<Real>("n_g", "gamma refractive index");
  params.addRequiredCoupledVar("var1", "the change in this refractive index");
  return params;
}


RefractiveIndex::RefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _na(getParam<Real>("n_a")),
   _nb(getParam<Real>("n_b")),
   _ng(getParam<Real>("n_g")),
  _var1(coupledValue("var1"))
{
}

Real
RefractiveIndex::computeValue()
{
  // the diagonals are related to the B1, B2, B3 terms in rotated indicatrix
//std::pow(  (1.0 / ( _beta_tensor[_qp](_index_i, _index_j)  ) ), 3.0)
  if (_component == 0)
    return _na - _var1[_qp];
  else if (_component == 1)
    return _nb - _var1[_qp];
  else if (_component == 2)
    return _ng - _var1[_qp];
  else
    return 0.0;
}


