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

#include "DielectricTensor.h"
#include "ComputeElectrostrictiveTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", DielectricTensor);

InputParameters DielectricTensor::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the local dielectric constant given by the classic (and likely wrong) thermodynamic relation.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0,  "The z component of the elastic displacement");
  params.addRequiredParam<Real>("alpha1", "_alpha1 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "_alpha11 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "_alpha12 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "_alpha111 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "_alpha112 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "_alpha123 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("first_deriv", "direction of first derivative");
  params.addRequiredParam<Real>("second_deriv", "direction of second derivative");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

DielectricTensor::DielectricTensor(const InputParameters & parameters) :
  AuxKernel(parameters),
  _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
  _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _disp_x_grad(coupledGradient("disp_x")),
  _disp_y_grad(coupledGradient("disp_y")),
  _disp_z_grad(coupledGradient("disp_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha12(getParam<Real>("alpha12")),
  _alpha111(getParam<Real>("alpha111")),
  _alpha112(getParam<Real>("alpha112")),
  _alpha123(getParam<Real>("alpha123")),
  _first_deriv(getParam<Real>("first_deriv")),
  _second_deriv(getParam<Real>("second_deriv")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
DielectricTensor::computeValue()
{
  //normal components:

  if (_first_deriv == 0 && _second_deriv == 0)
  {
    return 0.0;
  }
  else if (_first_deriv == 1 && _second_deriv == 1)
  {
    return 0.0;
  }
  else if (_first_deriv == 2 && _second_deriv == 2)
  {
    return 0.0;
  }
  //shears:
  else if (_first_deriv == 0 && _second_deriv == 1)
  {
    return 0.0;
  }
  else if (_first_deriv == 0 && _second_deriv == 2)
  {
    return 0.0;
  }
  else if (_first_deriv == 1 && _second_deriv == 2)
  {
    return 0.0;
  }
  else
    return 0.0;
}
