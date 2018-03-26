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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "RefractiveIndex.h"

registerMooseObject("FerretApp", RefractiveIndex);

template<>

InputParameters validParams<RefractiveIndex>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("index_j", "component");
  params.addRequiredParam<unsigned int>("index_k", "component");
  //We'll put some bools here to turn on and off different effects (elasto-, electro-, polar-,...)
  params.addParam<bool>("electro", false, "If this is true then electrooptic effect will be introduced");
  params.addParam<bool>("elasto", false, "If this is true then elastooptic effect will be introduced");
  params.addParam<bool>("polar", false, "If this is true then polaroptic effect will be introduced");
  params.addRequiredCoupledVar("var1", "the change in this refractive index");
  params.addCoupledVar("var2", 0.0, "the change in this refractive index");
  return params;
}


RefractiveIndex::RefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_j(getParam<unsigned int>("index_j")),
   _index_k(getParam<unsigned int>("index_k")),
   _electro(parameters.get<bool>("electro")),
   _elasto(parameters.get<bool>("elasto")),
   _polar(parameters.get<bool>("polar")),
   _indicatrix(getMaterialProperty<RankTwoTensor>("indicatrix")),
   _var1(coupledValue("var1")),
   _var2(coupledValue("var2"))
{
}

Real
RefractiveIndex::computeValue()
{
  // the diagonals are related to the B1, B2, B3 terms in rotated indicatrix
  if (_electro == true && _elasto == false && _polar == false)
    return std::pow(  (1.0 / ( _indicatrix[_qp](_index_j, _index_k)  ) ), 0.5) + _var1[_qp];
  else if (_elasto == true && _electro == false && _polar == false)
    return std::pow(  (1.0 / ( _indicatrix[_qp](_index_j, _index_k)  ) ), 0.5) + _var1[_qp];
  else if (_elasto == true && _electro == true && _polar == false)
    return std::pow(  (1.0 / ( _indicatrix[_qp](_index_j, _index_k)  ) ), 0.5) + _var1[_qp] + _var2[_qp];
  else
    return 0.0;
}
