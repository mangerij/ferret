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

#include "JacobiansRotoBulkEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", JacobiansRotoBulkEnergy);

InputParameters JacobiansRotoBulkEnergy::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the jacobian entries due to the microforce from the (roto) bulk terms.");
  params.addRequiredParam<unsigned int>("index_i", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<unsigned int>("index_j", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt");
  return params;
}

JacobiansRotoBulkEnergy::JacobiansRotoBulkEnergy(const InputParameters & parameters) :
  AuxKernel(parameters),
  _index_i(getParam<unsigned int>("index_i")),
  _index_j(getParam<unsigned int>("index_j")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
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
JacobiansRotoBulkEnergy::computeValue()
{
  if (_index_i == 0)
  {
    if (_index_j == 0) //on-diag
    {
      return (2*_beta1[_qp] + 12*_beta11[_qp]*Utility::pow<2>(_antiphase_A_x[_qp]) + 30*_beta111[_qp]*Utility::pow<4>(_antiphase_A_x[_qp]) + 56*_beta1111[_qp]*Utility::pow<6>(_antiphase_A_x[_qp]) + 
   2*_beta123[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + _beta12[_qp]*(2*Utility::pow<2>(_antiphase_A_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_z[_qp])) + 
   _beta1122[_qp]*(12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_antiphase_A_y[_qp]) + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_antiphase_A_z[_qp])) + _beta1123[_qp]*(12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + 2*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + 
      2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_antiphase_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiphase_A_y[_qp]) + 2*Utility::pow<4>(_antiphase_A_z[_qp]) + 
      12*Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiphase_A_y[_qp]) + 2*Utility::pow<6>(_antiphase_A_z[_qp]) + 
      30*Utility::pow<4>(_antiphase_A_x[_qp])*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))));
    }
    else if (_index_j == 1)
    {
      return (4*_beta12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_antiphase_A_y[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp] + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_y[_qp])) + 
   _beta1112[_qp]*(12*Utility::pow<5>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_antiphase_A_y[_qp])) + 4*_beta123[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp]) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp]) + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + 
      4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<4>(_antiphase_A_z[_qp])));
    }
    else if (_index_j == 2)
    {
      return (4*_beta12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_z[_qp] + 4*_beta123[_qp]*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_antiphase_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp] + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_z[_qp])) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<4>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 
      8*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<3>(_antiphase_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_antiphase_A_z[_qp])));
    }
    else
      return 0.0;
  }
  else if (_index_i == 1)
  {
    if (_index_j == 0)
    {
      return (4*_beta12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_antiphase_A_y[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp] + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_y[_qp])) + 
   _beta1112[_qp]*(12*Utility::pow<5>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_antiphase_A_y[_qp])) + 4*_beta123[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp]) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp]) + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + 
      4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<4>(_antiphase_A_z[_qp])));
    }
    else if (_index_j == 1) //on-diag
    {
      return (2*_beta1[_qp] + 12*_beta11[_qp]*Utility::pow<2>(_antiphase_A_y[_qp]) + 30*_beta111[_qp]*Utility::pow<4>(_antiphase_A_y[_qp]) + 56*_beta1111[_qp]*Utility::pow<6>(_antiphase_A_y[_qp]) + 
   2*_beta123[_qp]*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + _beta12[_qp]*(2*Utility::pow<2>(_antiphase_A_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_z[_qp])) + 
   _beta1123[_qp]*(2*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + 
      2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_antiphase_A_z[_qp])) + 
   _beta1122[_qp]*(12*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp]) + 12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_antiphase_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiphase_A_x[_qp]) + 2*Utility::pow<4>(_antiphase_A_z[_qp]) + 
      12*Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiphase_A_x[_qp]) + 2*Utility::pow<6>(_antiphase_A_z[_qp]) + 
      30*Utility::pow<4>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))));
    }
    else if (_index_j == 2)
    {
      return (4*_beta12[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp] + 4*_beta123[_qp]*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_antiphase_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 8*_antiphase_A_y[_qp]*Utility::pow<3>(_antiphase_A_z[_qp])) + 
   _beta1123[_qp]*(4*Utility::pow<4>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp] + 8*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<3>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 
      8*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*Utility::pow<3>(_antiphase_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_antiphase_A_z[_qp])));
    }
    else
      return 0.0;
  }
  else if (_index_i == 2)
  {
    if (_index_j == 0)
    {
      return (4*_beta12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_z[_qp] + 4*_beta123[_qp]*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_antiphase_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp] + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_z[_qp])) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<4>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 
      8*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<3>(_antiphase_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_antiphase_A_z[_qp])));
    }
    else if (_index_j == 1)
    {
      return (4*_beta12[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp] + 4*_beta123[_qp]*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_antiphase_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 8*_antiphase_A_y[_qp]*Utility::pow<3>(_antiphase_A_z[_qp])) + 
   _beta1123[_qp]*(4*Utility::pow<4>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp] + 8*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<3>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 
      8*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*Utility::pow<3>(_antiphase_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_antiphase_A_z[_qp])));
    }
    else if (_index_j == 2) //on-diag
    {
      return (2*_beta1[_qp] + 2*_beta123[_qp]*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp]) + _beta12[_qp]*(2*Utility::pow<2>(_antiphase_A_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])) + 
   12*_beta11[_qp]*Utility::pow<2>(_antiphase_A_z[_qp]) + 30*_beta111[_qp]*Utility::pow<4>(_antiphase_A_z[_qp]) + 56*_beta1111[_qp]*Utility::pow<6>(_antiphase_A_z[_qp]) + 
   _beta1123[_qp]*(2*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_antiphase_A_y[_qp]) + 
      12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])) + 
   _beta1122[_qp]*(12*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + 12*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiphase_A_x[_qp]) + 2*Utility::pow<4>(_antiphase_A_y[_qp]) + 
      12*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_antiphase_A_z[_qp])) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiphase_A_x[_qp]) + 2*Utility::pow<6>(_antiphase_A_y[_qp]) + 
      30*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<4>(_antiphase_A_z[_qp])));
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
