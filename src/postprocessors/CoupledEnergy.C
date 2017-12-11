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

#include "CoupledEnergy.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ComputeEigenstrain.h"
#include "libmesh/utility.h"


template<>
InputParameters validParams<CoupledEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("disp_x", "The x component of the elasticity displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elasticity displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elasticity displacement");
  params.addParam<Real>("artificial", 1.0, "term used to artificially increase coupling");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

CoupledEnergy::CoupledEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
  _eigenstrain(getMaterialProperty<RankTwoTensor>("eigenstrain")),
  _disp_x_grad(coupledGradient("disp_x")),
  _disp_y_grad(coupledGradient("disp_y")),
  _disp_z_grad(coupledGradient("disp_z")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _artificial(getParam<Real>("artificial")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
CoupledEnergy::computeQpIntegral()
{
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  Real sum3 = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  RealVectorValue v0(_eigenstrain[_qp](0,0), _eigenstrain[_qp](0,1), _eigenstrain[_qp](0,2));
  RealVectorValue v1(_eigenstrain[_qp](1,0), _eigenstrain[_qp](1,1), _eigenstrain[_qp](1,2));
  RealVectorValue v2(_eigenstrain[_qp](2,0), _eigenstrain[_qp](2,1), _eigenstrain[_qp](2,2));

  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], 0, w);
  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], 0, w);
  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], 0, w);

  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], 1, w);
  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], 1, w);
  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], 1, w);

  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], 2, w);
  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], 2, w);
  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], 2, w);

  return - _artificial * Utility::pow<3>(_len_scale) * ( sum1 * _polar_x[_qp] + sum2 * _polar_y[_qp] + sum3 * _polar_z[_qp]);
}
