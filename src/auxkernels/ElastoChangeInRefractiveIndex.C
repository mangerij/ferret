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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ElastoChangeInRefractiveIndex.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ElastoChangeInRefractiveIndex);

InputParameters ElastoChangeInRefractiveIndex::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the changes to local refractive index due to the elastooptic effect.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding changes in the refractive index's index (0, 1, 2, 3, 4, 5)");
  params.addRequiredCoupledVar("u_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the elastic displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the elastic displacement");
  return params;
}


ElastoChangeInRefractiveIndex::ElastoChangeInRefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _n1(getMaterialProperty<Real>("n1")),
   _n2(getMaterialProperty<Real>("n2")),
   _n3(getMaterialProperty<Real>("n3")),
   _n4(getMaterialProperty<Real>("n4")),
   _n5(getMaterialProperty<Real>("n5")),
   _n6(getMaterialProperty<Real>("n6")),
   _p1111(getMaterialProperty<Real>("p1111")),
   _p1122(getMaterialProperty<Real>("p1122")),
   _p1212(getMaterialProperty<Real>("p1212"))
{
}

Real
ElastoChangeInRefractiveIndex::computeValue()
{
  // consort with pp. 247 in Nye, Equation (19) and (28). The expressions for each of the six components of the refractive index change should be:
  // note that the impermeability tensor and the refractive index should be related by some matrix algebra...
  // Our focus is [001] oriented ferroelectrics _only_.
  //
  if (_component == 0)
  {
    return - 0.5 * Utility::pow<3>(_n1[_qp]) *  (_p1111[_qp]*_u_x_grad[_qp](0) + _p1122[_qp]*_u_y_grad[_qp](1) + _p1122[_qp]*_u_z_grad[_qp](2));
  }
  else if (_component == 1)
  {
    return - 0.5 * Utility::pow<3>(_n2[_qp]) *  (_p1122[_qp]*_u_x_grad[_qp](0) + _p1111[_qp]*_u_y_grad[_qp](1) + _p1122[_qp]*_u_z_grad[_qp](2));
  }
  else if (_component == 2)
  {
    return - 0.5 * Utility::pow<3>(_n3[_qp]) *  (_p1122[_qp]*_u_x_grad[_qp](0) + _p1122[_qp]*_u_y_grad[_qp](1) + _p1111[_qp]*_u_z_grad[_qp](2));
  }
  else if (_component == 3)
  {
    return - Utility::pow<3>(_n4[_qp]) *  (_p1212[_qp]*_u_y_grad[_qp](2));
  }
  else if (_component == 4)
  {
    return - Utility::pow<3>(_n5[_qp]) *  (_p1212[_qp]*_u_z_grad[_qp](0));
  }
  else if (_component == 5)
  {
    return - Utility::pow<3>(_n6[_qp]) *  (_p1212[_qp]*_u_x_grad[_qp](1));
  }
  else
    return 0.0;
}
