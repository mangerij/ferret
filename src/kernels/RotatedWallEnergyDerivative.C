/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "RotatedWallEnergyDerivative.h"

template<>
InputParameters validParams<RotatedWallEnergyDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  params.addRequiredParam<Real>("G110", "Domain wall coefficient");
  params.addRequiredParam<Real>("G11/G110", "Domain wall coefficient ratio");
  params.addRequiredParam<Real>("G12/G110", "Domain wall coefficient ratio");
  params.addRequiredParam<Real>("G44/G110", "Domain wall coefficient ratio");
  params.addRequiredParam<Real>("G44P/G110", "Domain wall coefficient ratio");
  ///params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

RotatedWallEnergyDerivative::RotatedWallEnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_x_var(coupled("polar_x")),
  _polar_y_var(coupled("polar_y")),
  _polar_z_var(coupled("polar_z")),
  _polar_i_grad((_component==0)? coupledGradient("polar_x") :(_component==1)? coupledGradient("polar_y"): coupledGradient("polar_z")),
  _polar_j_grad((_component==0)? coupledGradient("polar_y"): (_component==1)? coupledGradient("polar_z"): coupledGradient("polar_x")),
  _polar_k_grad((_component==0)? coupledGradient("polar_z"): (_component==1)? coupledGradient("polar_x"): coupledGradient("polar_y")),
  _Euler_angles(getParam<Real>("euler_angle_1"),
                getParam<Real>("euler_angle_2"),
                getParam<Real>("euler_angle_3")),
  _ii(_component),
  _jj((_component==0)? 1 : (_component==1)? 2: 0),
  _kk((_component==0)? 2 : (_component==1)? 0: 1),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11/G110") * _G110),
  _G12(getParam<Real>("G12/G110") * _G110),
  _G44(getParam<Real>("G44/G110") * _G110),
  _G44P(getParam<Real>("G44P/G110") * _G110),
  _len_scale(getParam<Real>("len_scale"))
{
  ///only for debug purpose
  std::cout<<"_G110 = "<<_G110<<"\n";
  std::cout<<"_G11 ="<<_G11<<"\n";
  std::cout<<"_G12 = "<<_G12<<"\n";
  std::cout<<"_G44 = "<<_G44<<"\n";
  std::cout<<"_G44P = "<<_G44P<<"\n";



}

Real
RotatedWallEnergyDerivative::computeQpResidual()
{
  /// need if statements here to construct the grad polar vector rotation:
  if(_component == 0)
    { //_polar_i_grad[_qp] = polar_x_grad, j = y, k = z
      Real Rwall = 0.0;
      RotationTensor R(_Euler_angles);
      RealVectorValue pgi = R(0,0) * _polar_i_grad[_qp] + R(0,1) * _polar_j_grad[_qp] + R(0,2) * _polar_k_grad[_qp]; 
      RealVectorValue pgj = R(1,0) * _polar_i_grad[_qp] + R(1,1) * _polar_j_grad[_qp] + R(1,2) * _polar_k_grad[_qp];
      RealVectorValue pgk = R(2,0) * _polar_i_grad[_qp] + R(2,1) * _polar_j_grad[_qp] + R(2,2) * _polar_k_grad[_qp];
      Rwall += (_G11 * pgi(_ii) * _grad_test[_i][_qp](_ii) + _G12 * (pgj(_jj) + pgk(_kk)) * _grad_test[_i][_qp](_ii) +
      _G44 * (pgi(_jj) + pgj(_ii)) * _grad_test[_i][_qp](_jj) + _G44 * (pgi(_kk)+pgk(_ii)) * _grad_test[_i][_qp](_kk) +
	   _G44P * (pgi(_jj) - pgj(_ii)) * _grad_test[_i][_qp](_jj) + _G44P * (pgi(_kk) - pgk(_ii)) * _grad_test[_i][_qp](_kk)) * _len_scale;
      ///  Moose::out << "\n R_wall-"; std::cout << _component << " = " << Rwall;
      return Rwall;
    }
  else if(_component == 1)
    { //_polar_i_grad[_qp] = polar_y_grad, j = z, k = x
      Real Rwall = 0.0;
      RotationTensor R(_Euler_angles);
      RealVectorValue pgi = R(0,0) * _polar_k_grad[_qp] + R(0,1) * _polar_i_grad[_qp] + R(0,2) * _polar_j_grad[_qp]; 
      RealVectorValue pgj = R(1,0) * _polar_k_grad[_qp] + R(1,1) * _polar_i_grad[_qp] + R(1,2) * _polar_j_grad[_qp];
      RealVectorValue pgk = R(2,0) * _polar_k_grad[_qp] + R(2,1) * _polar_i_grad[_qp] + R(2,2) * _polar_j_grad[_qp];
      Rwall += (_G11 * pgi(_ii) * _grad_test[_i][_qp](_ii) + _G12 * (pgj(_jj) + pgk(_kk)) * _grad_test[_i][_qp](_ii) +
      _G44 * (pgi(_jj) + pgj(_ii)) * _grad_test[_i][_qp](_jj) + _G44 * (pgi(_kk)+pgk(_ii)) * _grad_test[_i][_qp](_kk) +
	   _G44P * (pgi(_jj) - pgj(_ii)) * _grad_test[_i][_qp](_jj) + _G44P * (pgi(_kk) - pgk(_ii)) * _grad_test[_i][_qp](_kk)) * _len_scale;
      ///  Moose::out << "\n R_wall-"; std::cout << _component << " = " << Rwall;
      return Rwall;
    }
  else
    { //_polar_i_grad[_qp] = polar_z_grad, j = x, k = y
      Real Rwall = 0.0;
      RotationTensor R(_Euler_angles);
      RealVectorValue pgi = R(0,0) * _polar_j_grad[_qp] + R(0,1) * _polar_k_grad[_qp] + R(0,2) * _polar_i_grad[_qp]; 
      RealVectorValue pgj = R(1,0) * _polar_j_grad[_qp] + R(1,1) * _polar_k_grad[_qp] + R(1,2) * _polar_i_grad[_qp];
      RealVectorValue pgk = R(2,0) * _polar_j_grad[_qp] + R(2,1) * _polar_k_grad[_qp] + R(2,2) * _polar_i_grad[_qp];
      Rwall += (_G11 * pgi(_ii) * _grad_test[_i][_qp](_ii) + _G12 * (pgj(_jj) + pgk(_kk)) * _grad_test[_i][_qp](_ii) +
      _G44 * (pgi(_jj) + pgj(_ii)) * _grad_test[_i][_qp](_jj) + _G44 * (pgi(_kk) + pgk(_ii)) * _grad_test[_i][_qp](_kk) +
	   _G44P * (pgi(_jj) - pgj(_ii)) * _grad_test[_i][_qp](_jj) + _G44P * (pgi(_kk) - pgk(_ii)) * _grad_test[_i][_qp](_kk)) * _len_scale;
      ///  Moose::out << "\n R_wall-"; std::cout << _component << " = " << Rwall;
      return Rwall;
    }
}

Real
RotatedWallEnergyDerivative::computeQpJacobian()
{
  return (_G11 * _grad_phi[_j][_qp](_ii) * _grad_test[_i][_qp](_ii) +
    (_G44 + _G44P) * _grad_phi[_j][_qp](_jj) * _grad_test[_i][_qp](_jj) +
	  (_G44 + _G44P) * _grad_phi[_j][_qp](_kk) * _grad_test[_i][_qp](_kk)) * _len_scale;
}

Real
RotatedWallEnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  //mooseAssert(jvar!=variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var)
  {
    const unsigned int _jj = (jvar==_polar_x_var)? 0: (jvar==_polar_y_var)? 1 : 2;
    return (_G12 * _grad_phi[_j][_qp](_jj) * _grad_test[_i][_qp](_ii) + (_G44 - _G44P) * _grad_phi[_j][_qp](_ii) * _grad_test[_i][_qp](_jj)) * _len_scale;
  }
  else
  {
    return 0.0;
  }
}
