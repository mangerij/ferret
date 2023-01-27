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

#include "JacobiansRotopolarCoupledEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", JacobiansRotopolarCoupledEnergy);

InputParameters JacobiansRotopolarCoupledEnergy::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the jacobian entries for the rotopolar microforce.");
  params.addRequiredParam<unsigned int>("index_i", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<unsigned int>("index_j", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

JacobiansRotopolarCoupledEnergy::JacobiansRotopolarCoupledEnergy(const InputParameters & parameters) :
  AuxKernel(parameters),
  _index_i(getParam<unsigned int>("index_i")),
  _index_j(getParam<unsigned int>("index_j")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _antiphase_A_x(coupledValue("antiphase_A_x")),
  _antiphase_A_y(coupledValue("antiphase_A_y")),
  _antiphase_A_z(coupledValue("antiphase_A_z")),
  _t1111(getMaterialProperty<Real>("t1111")),
  _t1122(getMaterialProperty<Real>("t1122")),
  _t1212(getMaterialProperty<Real>("t1212")),
  _t42111111(getMaterialProperty<Real>("t42111111")),
  _t24111111(getMaterialProperty<Real>("t24111111")),
  _t42111122(getMaterialProperty<Real>("t42111122")),
  _t24112222(getMaterialProperty<Real>("t24112222")),
  _t42112233(getMaterialProperty<Real>("t42112233")),
  _t24112233(getMaterialProperty<Real>("t24112233")),
  _t42112211(getMaterialProperty<Real>("t42112211")),
  _t24111122(getMaterialProperty<Real>("t24111122")),
  _t42111212(getMaterialProperty<Real>("t42111212")),
  _t42123312(getMaterialProperty<Real>("t42123312")),
  _t24121112(getMaterialProperty<Real>("t24121112")),
  _t24121233(getMaterialProperty<Real>("t24121233")),
  _t6211111111(getMaterialProperty<Real>("t6211111111")),
  _t2611111111(getMaterialProperty<Real>("t2611111111")),
  _t6211111122(getMaterialProperty<Real>("t6211111122")),
  _t2611222222(getMaterialProperty<Real>("t2611222222")),
  _t4411111111(getMaterialProperty<Real>("t4411111111")),
  _t4411112222(getMaterialProperty<Real>("t4411112222"))
{
}

Real
JacobiansRotopolarCoupledEnergy::computeValue()
{
  if (_index_i == 0)
  {
    if (_index_j == 0) //on-diag
    {
      return (2*Utility::pow<2>(_antiphase_A_x[_qp])*_t1111[_qp] + 2*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_t1122[_qp] + 
   2*Utility::pow<4>(_antiphase_A_x[_qp])*_t24111111[_qp] + 2*Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_t24111122[_qp] + 
   2*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*_t24112222[_qp] + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])*_t24112233[_qp] + 
   2*Utility::pow<6>(_antiphase_A_x[_qp])*_t2611111111[_qp] + 2*(Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp]))*_t2611222222[_qp] + 
   12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t42111111[_qp] + 
   12*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp])*_t42111122[_qp] + 
   6*_polar_x[_qp]*(_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t42111212[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      2*Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*_t42112211[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42112233[_qp] + 
   2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp]*_t42123312[_qp] + 12*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t4411111111[_qp] + 
   12*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp])*_t4411112222[_qp] + 
   30*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp])*_t6211111111[_qp] + 
   30*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_x[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 1)
    {
      return (_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_t1212[_qp] + (Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_y[_qp]))*_t24121112[_qp] + 
   _antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_t24121233[_qp] + (3*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 3*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*
    _t42111212[_qp] + (4*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t42112211[_qp] + 
   4*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp]*_t42112233[_qp] + 
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp]);
    }
    else if (_index_j == 2)
    {
      return (_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_t1212[_qp] + (Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_z[_qp]))*_t24121112[_qp] + 
   _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_t24121233[_qp] + (3*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 3*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*
    _t42111212[_qp] + (4*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 4*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t42112211[_qp] + 
   4*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp]*_t42112233[_qp] + 
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (_index_j == 3)
    {
      return (4*_antiphase_A_x[_qp]*_polar_x[_qp]*_t1111[_qp] + (_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp]*_t24111111[_qp] + 
   4*_antiphase_A_x[_qp]*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_polar_x[_qp]*_t24111122[_qp] + 
   (Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_x[_qp])*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_x[_qp]*_t2611111111[_qp] + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111111[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + _antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_x[_qp])*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_x[_qp]*_polar_x[_qp]*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] + 
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411111111[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111111[_qp]);
    }
    else if (_index_j == 4)
    {
      return (4*_antiphase_A_y[_qp]*_polar_x[_qp]*_t1122[_qp] + _antiphase_A_x[_qp]*_polar_y[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_x[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_x[_qp]*_t24112222[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t24112233[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_y[_qp] + 3*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp])*_t24121112[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_x[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122[_qp] + (3*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112211[_qp] + 4*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 
   (2*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 5)
    {
      return (4*_antiphase_A_z[_qp]*_polar_x[_qp]*_t1122[_qp] + _antiphase_A_x[_qp]*_polar_z[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_t24112233[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_z[_qp] + 3*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122[_qp] + (3*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] + 4*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112233[_qp] + 
   (2*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122[_qp]);
    }
    else
      return 0.0;
  }
  else if (_index_i == 1)
  {
    if (_index_j == 0)
    {
      return (_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_t1212[_qp] + (Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_y[_qp]))*_t24121112[_qp] + 
   _antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_t24121233[_qp] + (3*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 3*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*
    _t42111212[_qp] + (4*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t42112211[_qp] + 
   4*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp]*_t42112233[_qp] + 
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp]);
    }
    else if (_index_j == 1) //on-diag
    {
      return (2*Utility::pow<2>(_antiphase_A_y[_qp])*_t1111[_qp] + 2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_t1122[_qp] + 
   2*Utility::pow<4>(_antiphase_A_y[_qp])*_t24111111[_qp] + 2*Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_t24111122[_qp] + 
   2*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*_t24112222[_qp] + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])*_t24112233[_qp] + 
   2*Utility::pow<6>(_antiphase_A_y[_qp])*_t2611111111[_qp] + 2*(Utility::pow<6>(_antiphase_A_x[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp]))*_t2611222222[_qp] + 
   12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t42111111[_qp] + 
   12*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp])*_t42111122[_qp] + 
   6*_polar_y[_qp]*(_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t42111212[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      2*Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])))*_t42112211[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42112233[_qp] + 
   2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp]*_t42123312[_qp] + 12*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t4411111111[_qp] + 
   12*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp])*_t4411112222[_qp] + 
   30*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp])*_t6211111111[_qp] + 
   30*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_y[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 2)
    {
      return (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_t1212[_qp] + (Utility::pow<3>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + _antiphase_A_y[_qp]*Utility::pow<3>(_antiphase_A_z[_qp]))*_t24121112[_qp] + 
   Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_t24121233[_qp] + (3*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 3*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*
    _t42111212[_qp] + (4*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 4*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t42112211[_qp] + 
   4*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42112233[_qp] + 
   (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (_index_j == 3)
    {
      return (4*_antiphase_A_x[_qp]*_polar_y[_qp]*_t1122[_qp] + _antiphase_A_y[_qp]*_polar_x[_qp]*_t1212[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_y[_qp]*_t24112222[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_x[_qp])*_t24121112[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_y[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122[_qp] + (_antiphase_A_y[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112211[_qp] + 4*_antiphase_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 
   (2*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411112222[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 4)
    {
      return (4*_antiphase_A_y[_qp]*_polar_y[_qp]*_t1111[_qp] + (_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t24111111[_qp] + 
   4*_antiphase_A_y[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_polar_y[_qp]*_t24111122[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t2611111111[_qp] + 8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111111[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_y[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_y[_qp]*_polar_y[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] + 
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411111111[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111111[_qp]);
    }
    else if (_index_j == 5)
    {
      return (4*_antiphase_A_z[_qp]*_polar_y[_qp]*_t1122[_qp] + _antiphase_A_y[_qp]*_polar_z[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp]*_t24112233[_qp] + 
   (Utility::pow<3>(_antiphase_A_y[_qp])*_polar_z[_qp] + 3*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122[_qp] + (3*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] + 4*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112233[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411112222[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122[_qp]);
    }
    else
      return 0.0;
  }
  else if (_index_i == 2)
  {
    if (_index_j == 0)
    {
      return (_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_t1212[_qp] + (Utility::pow<3>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_antiphase_A_z[_qp]))*_t24121112[_qp] + 
   _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_t24121233[_qp] + (3*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 3*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*
    _t42111212[_qp] + (4*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 4*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t42112211[_qp] + 
   4*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp]*_t42112233[_qp] + 
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (_index_j == 1)
    {
      return (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_t1212[_qp] + (Utility::pow<3>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp] + _antiphase_A_y[_qp]*Utility::pow<3>(_antiphase_A_z[_qp]))*_t24121112[_qp] + 
   Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_t24121233[_qp] + (3*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 3*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*
    _t42111212[_qp] + (4*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 4*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t42112211[_qp] + 
   4*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42112233[_qp] + 
   (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (_index_j == 2)  //on-diag
    {
      return (2*Utility::pow<2>(_antiphase_A_z[_qp])*_t1111[_qp] + 2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_t1122[_qp] + 
   2*Utility::pow<4>(_antiphase_A_z[_qp])*_t24111111[_qp] + 2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_antiphase_A_z[_qp])*_t24111122[_qp] + 
   2*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]))*_t24112222[_qp] + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*_t24112233[_qp] + 
   2*Utility::pow<6>(_antiphase_A_z[_qp])*_t2611111111[_qp] + 2*(Utility::pow<6>(_antiphase_A_x[_qp]) + Utility::pow<6>(_antiphase_A_y[_qp]))*_t2611222222[_qp] + 
   12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42111111[_qp] + 
   12*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp])*_t42111122[_qp] + 
   6*(_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp])*_polar_z[_qp]*_t42111212[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 
      2*Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp])))*_t42112211[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t42112233[_qp] + 
   2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_t42123312[_qp] + 12*Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t4411111111[_qp] + 
   12*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp])*_t4411112222[_qp] + 
   30*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp])*_t6211111111[_qp] + 
   30*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<4>(_polar_z[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 3)
    {
      return (4*_antiphase_A_x[_qp]*_polar_z[_qp]*_t1122[_qp] + _antiphase_A_z[_qp]*_polar_x[_qp]*_t1212[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_z[_qp]*_t24112222[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_x[_qp])*_t24121112[_qp] + 
   (Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_z[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122[_qp] + (_antiphase_A_z[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112211[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112233[_qp] + 
   (_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 4)
    {
      return (4*_antiphase_A_y[_qp]*_polar_z[_qp]*_t1122[_qp] + _antiphase_A_z[_qp]*_polar_y[_qp]*_t1212[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_z[_qp]*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_y[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122[_qp] + (_antiphase_A_z[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*_antiphase_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112211[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112233[_qp] + 
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 5)
    {
      return (4*_antiphase_A_z[_qp]*_polar_z[_qp]*_t1111[_qp] + (_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111111[_qp] + 
   4*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_antiphase_A_z[_qp]*_polar_z[_qp]*_t24111122[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp]))*_t24121112[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_y[_qp])*_t24121233[_qp] + 
   12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t2611111111[_qp] + 8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111111[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_z[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*_polar_z[_qp]*_t42112211[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411111111[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111111[_qp]);
    }
    else
      return 0.0;
  }
  else if (_index_i == 3)
  {
    if (_index_j == 0)
    {
      return (4*_antiphase_A_x[_qp]*_polar_x[_qp]*_t1111[_qp] + (_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp]*_t24111111[_qp] + 
   4*_antiphase_A_x[_qp]*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_polar_x[_qp]*_t24111122[_qp] + 
   (Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_x[_qp])*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_x[_qp]*_t2611111111[_qp] + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111111[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + _antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_x[_qp])*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_x[_qp]*_polar_x[_qp]*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] + 
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411111111[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111111[_qp]);
    }
    else if (_index_j == 1)
    {
      return (4*_antiphase_A_x[_qp]*_polar_y[_qp]*_t1122[_qp] + _antiphase_A_y[_qp]*_polar_x[_qp]*_t1212[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_y[_qp]*_t24112222[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_x[_qp])*_t24121112[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_y[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122[_qp] + (_antiphase_A_y[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112211[_qp] + 4*_antiphase_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 
   (2*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411112222[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 2)
    {
      return (4*_antiphase_A_x[_qp]*_polar_z[_qp]*_t1122[_qp] + _antiphase_A_z[_qp]*_polar_x[_qp]*_t1212[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_z[_qp]*_t24112222[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_x[_qp])*_t24121112[_qp] + 
   (Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_z[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122[_qp] + (_antiphase_A_z[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112211[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112233[_qp] + 
   (_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 3)  //on-diag
    {
      return (2*Utility::pow<2>(_polar_x[_qp])*_t1111[_qp] + (2*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))*_t1122[_qp] + 
   12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t24111111[_qp] + 
   (2*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 
      2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   (12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233[_qp] + 
   6*_antiphase_A_x[_qp]*_polar_x[_qp]*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t24121112[_qp] + 2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp]*_t24121233[_qp] + 
   30*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t2611111111[_qp] + 
   (30*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 30*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222[_qp] + 
   2*Utility::pow<4>(_polar_x[_qp])*_t42111111[_qp] + (2*Utility::pow<4>(_polar_y[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]))*_t42111122[_qp] + 
   2*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] + 
   2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp])*_t4411111111[_qp] + 
   (12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_y[_qp]) + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222[_qp] + 
   2*Utility::pow<6>(_polar_x[_qp])*_t6211111111[_qp] + (2*Utility::pow<6>(_polar_y[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]))*_t6211111122[_qp]);
    }
    else if (_index_j == 4)
    {
      return (_polar_x[_qp]*_polar_y[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24111122[_qp] + 
   4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t24121112[_qp] + 
   (Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   (Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212[_qp] + _polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (_index_j == 5)
    {
      return (_polar_x[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   (Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + _polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42123312[_qp]);
    }
    else
      return 0.0;
  }
  else if (_index_i == 4)
  {
    if (_index_j == 0)
    {
      return (4*_antiphase_A_y[_qp]*_polar_x[_qp]*_t1122[_qp] + _antiphase_A_x[_qp]*_polar_y[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_x[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_x[_qp]*_t24112222[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t24112233[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_y[_qp] + 3*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp])*_t24121112[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_x[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122[_qp] + (3*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112211[_qp] + 4*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 
   (2*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 1)
    {
      return (4*_antiphase_A_y[_qp]*_polar_y[_qp]*_t1111[_qp] + (_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t24111111[_qp] + 
   4*_antiphase_A_y[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_polar_y[_qp]*_t24111122[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t2611111111[_qp] + 8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111111[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_y[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_y[_qp]*_polar_y[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] + 
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411111111[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111111[_qp]);
    }
    else if (_index_j == 2)
    {
      return (4*_antiphase_A_y[_qp]*_polar_z[_qp]*_t1122[_qp] + _antiphase_A_z[_qp]*_polar_y[_qp]*_t1212[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_z[_qp]*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_y[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122[_qp] + (_antiphase_A_z[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*_antiphase_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112211[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112233[_qp] + 
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 3)
    {
      return (_polar_x[_qp]*_polar_y[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24111122[_qp] + 
   4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t24121112[_qp] + 
   (Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   (Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212[_qp] + _polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (_index_j == 4)  //on-diag
    {
      return (2*Utility::pow<2>(_polar_y[_qp])*_t1111[_qp] + (2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))*_t1122[_qp] + 
   12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t24111111[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + 
      2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   (12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233[_qp] + 
   6*_antiphase_A_y[_qp]*_polar_y[_qp]*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t24121112[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp]*_t24121233[_qp] + 
   30*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t2611111111[_qp] + 
   (30*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 30*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222[_qp] + 
   2*Utility::pow<4>(_polar_y[_qp])*_t42111111[_qp] + (2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]))*_t42111122[_qp] + 
   2*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] + 
   2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp])*_t4411111111[_qp] + 
   (12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222[_qp] + 
   2*Utility::pow<6>(_polar_y[_qp])*_t6211111111[_qp] + (2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]))*_t6211111122[_qp]);
    }
    else if (_index_j == 5)
    {
      return (_polar_y[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   (Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + _polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42123312[_qp]);
    }
    else
      return 0.0;
  }
  else if (_index_i == 5)
  {
    if (_index_j == 0)
    {
      return (4*_antiphase_A_z[_qp]*_polar_x[_qp]*_t1122[_qp] + _antiphase_A_x[_qp]*_polar_z[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_t24112233[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_z[_qp] + 3*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122[_qp] + (3*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] + 4*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112233[_qp] + 
   (2*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 1)
    {
      return (4*_antiphase_A_z[_qp]*_polar_y[_qp]*_t1122[_qp] + _antiphase_A_y[_qp]*_polar_z[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp]*_t24111122[_qp] + 
   8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp]*_t24112233[_qp] + 
   (Utility::pow<3>(_antiphase_A_y[_qp])*_polar_z[_qp] + 3*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t2611222222[_qp] + 
   8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122[_qp] + (3*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] + 4*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112233[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411112222[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122[_qp]);
    }
    else if (_index_j == 2)
    {
      return (4*_antiphase_A_z[_qp]*_polar_z[_qp]*_t1111[_qp] + (_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111111[_qp] + 
   4*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_antiphase_A_z[_qp]*_polar_z[_qp]*_t24111122[_qp] + 
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp]))*_t24121112[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_y[_qp])*_t24121233[_qp] + 
   12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t2611111111[_qp] + 8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111111[_qp] + 
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] + 
   4*_antiphase_A_z[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*_polar_z[_qp]*_t42112211[_qp] + 
   (_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42123312[_qp] + 
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411111111[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111111[_qp]);
    }
    else if (_index_j == 3)
    {
      return (_polar_x[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   (Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + _polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42123312[_qp]);
    }
    else if (_index_j == 4)
    {
      return (_polar_y[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_t24112233[_qp] + 
   (3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121112[_qp] + 
   (2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   (Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + _polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42123312[_qp]);
    }
    else if (_index_j == 5)  //on-diag
    {
      return (2*Utility::pow<2>(_polar_z[_qp])*_t1111[_qp] + (2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]))*_t1122[_qp] + 
   12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t24111111[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 
      2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   (12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112222[_qp] + 
   (2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112233[_qp] + 
   6*_antiphase_A_z[_qp]*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_polar_z[_qp]*_t24121112[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_t24121233[_qp] + 
   30*Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t2611111111[_qp] + 
   (30*Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 30*Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t2611222222[_qp] + 
   2*Utility::pow<4>(_polar_z[_qp])*_t42111111[_qp] + (2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_y[_qp]))*_t42111122[_qp] + 
   2*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] + 
   2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_t42112233[_qp] + 12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp])*_t4411111111[_qp] + 
   (12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_y[_qp]))*_t4411112222[_qp] + 
   2*Utility::pow<6>(_polar_z[_qp])*_t6211111111[_qp] + (2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_y[_qp]))*_t6211111122[_qp]);
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
