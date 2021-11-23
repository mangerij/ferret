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

#include "ReworkedRefractiveIndex.h"

registerMooseObject("FerretApp", ReworkedRefractiveIndex);

InputParameters ReworkedRefractiveIndex::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the refractive index's Voight index (0, 1, 2, 3, 4, 5)");
  params.addParam<bool>("electro", false, "If this is true then electrooptic effect will be introduced");
  params.addParam<bool>("elasto", false, "If this is true then elastooptic effect will be introduced");
  params.addParam<bool>("polar", false, "If this is true then polaroptic effect will be introduced");
  params.addRequiredCoupledVar("var1", "the change in this refractive index");
  params.addCoupledVar("var2", 0.0, "the change in this refractive index due to the second effect");
  params.addCoupledVar("var3", 0.0, "the change in this refractive index due to the third effect");
  return params;
}


ReworkedRefractiveIndex::ReworkedRefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _n1(getMaterialProperty<Real>("n1")),
   _n2(getMaterialProperty<Real>("n2")),
   _n3(getMaterialProperty<Real>("n3")),
   _n4(getMaterialProperty<Real>("n4")),
   _n5(getMaterialProperty<Real>("n5")),
   _n6(getMaterialProperty<Real>("n6")),
   _electro(parameters.get<bool>("electro")),
   _elasto(parameters.get<bool>("elasto")),
   _polar(parameters.get<bool>("polar")),
   _var1(coupledValue("var1")),
   _var2(coupledValue("var2")),
   _var3(coupledValue("var3"))
{
}

Real
ReworkedRefractiveIndex::computeValue()
{
  if (_component == 0)
  {
    if (_electro == true && _elasto == false && _polar == false)
      return _n1[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == false && _polar == false)
      return _n1[_qp] + _var1[_qp];
    else if (_elasto == false && _electro == false && _polar == true)
      return _n1[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == true && _polar == false)
      return _n1[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == false && _polar == true)
      return _n1[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == false && _electro == true && _polar == true)
      return _n1[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == true && _polar == true)
      return _n1[_qp] + _var1[_qp] + _var2[_qp] + _var3[_qp];
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (_electro == true && _elasto == false && _polar == false)
      return _n2[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == false && _polar == false)
      return _n2[_qp] + _var1[_qp];
    else if (_elasto == false && _electro == false && _polar == true)
      return _n2[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == true && _polar == false)
      return _n2[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == false && _polar == true)
      return _n2[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == false && _electro == true && _polar == true)
      return _n2[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == true && _polar == true)
      return _n2[_qp] + _var1[_qp] + _var2[_qp] + _var3[_qp];
    else
      return 0.0;
  }
  else if (_component == 2)
  {
    if (_electro == true && _elasto == false && _polar == false)
      return _n3[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == false && _polar == false)
      return _n3[_qp] + _var1[_qp];
    else if (_elasto == false && _electro == false && _polar == true)
      return _n3[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == true && _polar == false)
      return _n3[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == false && _polar == true)
      return _n3[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == false && _electro == true && _polar == true)
      return _n3[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == true && _polar == true)
      return _n3[_qp] + _var1[_qp] + _var2[_qp] + _var3[_qp];
    else
      return 0.0;
  }
  else if (_component == 3)
  {
    if (_electro == true && _elasto == false && _polar == false)
      return _n4[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == false && _polar == false)
      return _n4[_qp] + _var1[_qp];
    else if (_elasto == false && _electro == false && _polar == true)
      return _n4[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == true && _polar == false)
      return _n4[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == false && _polar == true)
      return _n4[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == false && _electro == true && _polar == true)
      return _n4[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == true && _polar == true)
      return _n4[_qp] + _var1[_qp] + _var2[_qp] + _var3[_qp];
    else
      return 0.0;
  }
  else if (_component == 4)
  {
    if (_electro == true && _elasto == false && _polar == false)
      return _n5[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == false && _polar == false)
      return _n5[_qp] + _var1[_qp];
    else if (_elasto == false && _electro == false && _polar == true)
      return _n5[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == true && _polar == false)
      return _n5[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == false && _polar == true)
      return _n5[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == false && _electro == true && _polar == true)
      return _n5[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == true && _polar == true)
      return _n5[_qp] + _var1[_qp] + _var2[_qp] + _var3[_qp];
    else
      return 0.0;
  }
  else if (_component == 5)
  {
    if (_electro == true && _elasto == false && _polar == false)
      return _n6[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == false && _polar == false)
      return _n6[_qp] + _var1[_qp];
    else if (_elasto == false && _electro == false && _polar == true)
      return _n6[_qp] + _var1[_qp];
    else if (_elasto == true && _electro == true && _polar == false)
      return _n6[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == false && _polar == true)
      return _n6[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == false && _electro == true && _polar == true)
      return _n6[_qp] + _var1[_qp] + _var2[_qp];
    else if (_elasto == true && _electro == true && _polar == true)
      return _n6[_qp] + _var1[_qp] + _var2[_qp] + _var3[_qp];
    else
      return 0.0;
  }
  else
    return 0.0;
}
