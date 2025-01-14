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

#include "DemagFieldAuxPML.h"

registerMooseObject("FerretApp", DemagFieldAuxPML);

InputParameters DemagFieldAuxPML::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Converts magnetostatic potential to the vector demagnetization field.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("deltasyminus", 1.e03, "scaled thickness of pml region");
  params.addParam<Real>("deltapyminus",1., "pole distance");
  params.addParam<Real>("deltawyminus",1.,"distance to be scaled (dimensionfull");
  params.addParam<Real>("y0pmlminus",-1., "position at which scaling will start");
  params.addCoupledVar("phi1", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", "The external magnetic potential variable (disabled for now)");
  return params;
}


DemagFieldAuxPML::DemagFieldAuxPML(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _deltasyminus(getParam<Real>("deltasyminus")),
   _deltapyminus(getParam<Real>("deltapyminus")),
   _deltawyminus(getParam<Real>("deltawyminus")),
  _y0pmlminus(getParam<Real>("y0pmlminus")),
   _phi1_grad(coupledGradient("phi1")),
   _potential_H_ext_grad(coupledGradient("potential_H_ext"))
{
}

Real
DemagFieldAuxPML::computeValue()

{
    if (_component == 1)
  {
    	 if  (_q_point[_qp](1) < _y0pmlminus)
	 {
	   const Real gamma = (_deltasyminus + _deltapyminus)/_deltasyminus;
	   const Real xi = -(_q_point[_qp](1)-_y0pmlminus)/_deltawyminus;
//	   const Real dudx = 1.+_deltapyminus/_deltawyminus*(xi/(gamma-xi))*((2.*gamma-xi)/(gamma-xi));
           const Real dudx = 1.+_deltapyminus/_deltawyminus*(1./(gamma-xi)*(1.+xi/(gamma-xi))-1./gamma);
           return - _phi1_grad[_qp](_component)/dudx - _potential_H_ext_grad[_qp](_component);
	  }
         else
	   return 0.0;
  }
    else
    return 0.0;
}
