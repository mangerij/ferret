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

#include "RotoBulkEnergyEighth.h"

registerMooseObject("FerretApp", RotoBulkEnergyEighth);

template<>
InputParameters validParams<RotoBulkEnergyEighth>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral whose integrand is the eighth order expansion of the AFD fields");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the AFD vector field");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the AFD vector field");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the AFD vector field");
  return params;
}

RotoBulkEnergyEighth::RotoBulkEnergyEighth(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _beta1(getMaterialProperty<Real>("beta1")),
   _beta11(getMaterialProperty<Real>("beta11")),
   _beta12(getMaterialProperty<Real>("beta12")),
   _beta111(getMaterialProperty<Real>("beta111")),
   _beta112(getMaterialProperty<Real>("beta112")),
   _beta123(getMaterialProperty<Real>("beta123")),
   _beta1111(getMaterialProperty<Real>("beta1111")),
   _beta1112(getMaterialProperty<Real>("beta1112")),
   _beta1122(getMaterialProperty<Real>("beta1122")),
   _beta1123(getMaterialProperty<Real>("beta1123"))
{
}

Real
RotoBulkEnergyEighth::computeQpIntegral()
{
  return (_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta1[_qp]*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp])) + _beta12[_qp]*(Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta11[_qp]*(Utility::pow<4>(_antiferrodis_A_x[_qp]) + Utility::pow<4>(_antiferrodis_A_y[_qp]) + Utility::pow<4>(_antiferrodis_A_z[_qp])) + _beta1123[_qp]*(Utility::pow<6>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1122[_qp]*(Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp]) + Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp]) + Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + _beta111[_qp]*(Utility::pow<6>(_antiferrodis_A_x[_qp]) + Utility::pow<6>(_antiferrodis_A_y[_qp]) + Utility::pow<6>(_antiferrodis_A_z[_qp])) + _beta1111[_qp]*(Utility::pow<8>(_antiferrodis_A_x[_qp]) + Utility::pow<8>(_antiferrodis_A_y[_qp]) + Utility::pow<8>(_antiferrodis_A_z[_qp])) + 
   _beta112[_qp]*((Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<4>(_antiferrodis_A_z[_qp]) + Utility::pow<4>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp])) + Utility::pow<4>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112[_qp]*((Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<6>(_antiferrodis_A_z[_qp]) + Utility::pow<6>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp])) + Utility::pow<6>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
}
