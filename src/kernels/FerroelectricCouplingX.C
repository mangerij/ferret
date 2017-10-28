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

#include "FerroelectricCouplingX.h"

class FerroelectricCouplingX;

template<>
InputParameters validParams<FerroelectricCouplingX>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("disp_x", "The x component of the displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("artificial", 1.0, "artificial increase coupling");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

FerroelectricCouplingX::FerroelectricCouplingX(const InputParameters & parameters)
  :Kernel(parameters),
   _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
   _component(getParam<unsigned int>("component")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FerroelectricCouplingX::computeQpResidual()
{
  Real sum = 0.0;
  Real Rp = 0.0;

  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], 0, p) * _polar_x[_qp];
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], 1, p) * _polar_y[_qp];
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], 2, p) * _polar_z[_qp];
  Rp = std::pow(_len_scale, 2.0) * sum;
  /// Moose::out << "\n R ="; std::cout << Rp;
  return Rp;
}


Real
FerroelectricCouplingX::computeQpJacobian()
{
  return 0.0; ///dRdu_i = 0 for this term!
}

Real
FerroelectricCouplingX::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component;
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  /// return = 0.0;
  if (jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
  {
    if (jvar == _polar_x_var)
    {
      coupled_component = 0;
      sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], coupled_component, p);
    }
    else if (jvar == _polar_y_var)
    {
      coupled_component = 1;
      sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], coupled_component, p);
    }
    else if (jvar == _polar_z_var)
    {
      coupled_component = 2;
      sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], coupled_component, p) ;
    }
    return std::pow(_len_scale, 2.0) * _phi[_j][_qp] * sum1;
  }
  else
  {
    return 0.0;
  }
}
