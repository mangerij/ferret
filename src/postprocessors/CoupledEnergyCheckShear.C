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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "CoupledEnergyCheckShear.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ComputeEigenstrain.h"

template<>
InputParameters validParams<CoupledEnergyCheckShear>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("disp_x", "The x component of the elasticity displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elasticity displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elasticity displacement");
  params.addParam<Real>("artificial", 1.0, "artificial increase coupling");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

CoupledEnergyCheckShear::CoupledEnergyCheckShear(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
  _stress_free_strain(getMaterialProperty<RankTwoTensor>("stress_free_strain")),
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
CoupledEnergyCheckShear::computeQpIntegral()
{
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  Real sum3 = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  RealVectorValue v0(_stress_free_strain[_qp](0,0), _stress_free_strain[_qp](0,1), _stress_free_strain[_qp](0,2));
  RealVectorValue v1(_stress_free_strain[_qp](1,0), _stress_free_strain[_qp](1,1), _stress_free_strain[_qp](1,2));
  RealVectorValue v2(_stress_free_strain[_qp](2,0), _stress_free_strain[_qp](2,1), _stress_free_strain[_qp](2,2));

  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, 0, w);
  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, 0, w);
  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, 0, w);

  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, 1, w);
  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, 1, w);
  sum2 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, 1, w);

  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, 2, w);
  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, 2, w);
  sum3 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, 2, w);

  return _artificial * std::pow(_len_scale, 3.0) * ( sum1 * _polar_x[_qp] + sum2 * _polar_y[_qp] + sum3 * _polar_z[_qp]) - 4 * _artificial * std::pow(_len_scale, 3.0) * (
         _electrostrictive_tensor[_qp](0,2,0,0) * w(0) * w(0) + _electrostrictive_tensor[_qp](0,2,1,0) * w(1) * w(0) + _electrostrictive_tensor[_qp](0,2,2,0) * w(2) * w(0) +
         _electrostrictive_tensor[_qp](0,2,0,1) * w(0) * w(1) + _electrostrictive_tensor[_qp](0,2,1,1) * w(1) * w(1) + _electrostrictive_tensor[_qp](0,2,2,1) * w(2) * w(1) +
         _electrostrictive_tensor[_qp](0,2,0,2) * w(0) * w(2) + _electrostrictive_tensor[_qp](0,2,1,2) * w(1) * w(2) + _electrostrictive_tensor[_qp](0,2,2,2) * w(2) * w(2)  );
}
