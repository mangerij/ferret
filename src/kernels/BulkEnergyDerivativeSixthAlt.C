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

#include "BulkEnergyDerivativeSixthAlt.h"
#include<cmath>

template<>
InputParameters validParams<BulkEnergyDerivativeSixthAlt>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

BulkEnergyDerivativeSixthAlt::BulkEnergyDerivativeSixthAlt(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _alpha1(getParam<Real>("alpha1")),
   _alpha11(getParam<Real>("alpha11")),
   _alpha12(getParam<Real>("alpha12")),
   _alpha111(getParam<Real>("alpha111")),
   _alpha112(getParam<Real>("alpha112")),
   _alpha123(getParam<Real>("alpha123")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkEnergyDerivativeSixthAlt::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_alpha1*_polar_x[_qp] + 4*_alpha11*std::pow(_polar_x[_qp],3) + 6*_alpha111*std::pow(_polar_x[_qp],5) + 
   2*_alpha123*_polar_x[_qp]*std::pow(_polar_y[_qp],2)*std::pow(_polar_z[_qp],2) + 
   _alpha12*(2*_polar_x[_qp]*std::pow(_polar_y[_qp],2) + 2*_polar_x[_qp]*std::pow(_polar_z[_qp],2)) + 
   _alpha112*(2*_polar_x[_qp]*std::pow(_polar_y[_qp],4) + 2*_polar_x[_qp]*std::pow(_polar_z[_qp],4) + 
      4*std::pow(_polar_x[_qp],3)*(std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*_alpha1*_polar_y[_qp] + 4*_alpha11*std::pow(_polar_y[_qp],3) + 6*_alpha111*std::pow(_polar_y[_qp],5) + 
   2*_alpha123*std::pow(_polar_x[_qp],2)*_polar_y[_qp]*std::pow(_polar_z[_qp],2) + 
   _alpha12*(2*std::pow(_polar_x[_qp],2)*_polar_y[_qp] + 2*_polar_y[_qp]*std::pow(_polar_z[_qp],2)) + 
   _alpha112*(2*std::pow(_polar_x[_qp],4)*_polar_y[_qp] + 2*_polar_y[_qp]*std::pow(_polar_z[_qp],4) + 
      4*std::pow(_polar_y[_qp],3)*(std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_alpha1*_polar_z[_qp] + 2*_alpha123*std::pow(_polar_x[_qp],2)*std::pow(_polar_y[_qp],2)*_polar_z[_qp] + 
   4*_alpha11*std::pow(_polar_z[_qp],3) + 6*_alpha111*std::pow(_polar_z[_qp],5) + 
   _alpha12*(2*std::pow(_polar_x[_qp],2)*_polar_z[_qp] + 2*std::pow(_polar_y[_qp],2)*_polar_z[_qp]) + 
   _alpha112*(2*std::pow(_polar_x[_qp],4)*_polar_z[_qp] + 2*std::pow(_polar_y[_qp],4)*_polar_z[_qp] + 
      4*(std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*std::pow(_polar_z[_qp],3)));
  }
  else
    return 0.0;
}

Real
BulkEnergyDerivativeSixthAlt::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_alpha1 + 12*_alpha11*std::pow(_polar_x[_qp],2) + 30*_alpha111*std::pow(_polar_x[_qp],4) + 
   2*_alpha123*std::pow(_polar_y[_qp],2)*std::pow(_polar_z[_qp],2) + 
   _alpha12*(2*std::pow(_polar_y[_qp],2) + 2*std::pow(_polar_z[_qp],2)) + 
   _alpha112*(2*std::pow(_polar_y[_qp],4) + 2*std::pow(_polar_z[_qp],4) + 
      12*std::pow(_polar_x[_qp],2)*(std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_alpha1 + 12*_alpha11*std::pow(_polar_y[_qp],2) + 30*_alpha111*std::pow(_polar_y[_qp],4) + 
   2*_alpha123*std::pow(_polar_x[_qp],2)*std::pow(_polar_z[_qp],2) + 
   _alpha12*(2*std::pow(_polar_x[_qp],2) + 2*std::pow(_polar_z[_qp],2)) + 
   _alpha112*(2*std::pow(_polar_x[_qp],4) + 2*std::pow(_polar_z[_qp],4) + 
      12*std::pow(_polar_y[_qp],2)*(std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_alpha1 + 2*_alpha123*std::pow(_polar_x[_qp],2)*std::pow(_polar_y[_qp],2) + 
   _alpha12*(2*std::pow(_polar_x[_qp],2) + 2*std::pow(_polar_y[_qp],2)) + 12*_alpha11*std::pow(_polar_z[_qp],2) + 
   30*_alpha111*std::pow(_polar_z[_qp],4) + _alpha112*
    (2*std::pow(_polar_x[_qp],4) + 2*std::pow(_polar_y[_qp],4) + 
      12*(std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*std::pow(_polar_z[_qp],2)));
  }
  else
    return 0.0;
}

Real
BulkEnergyDerivativeSixthAlt::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12*_polar_x[_qp]*_polar_y[_qp] + _alpha112*
    (8*std::pow(_polar_x[_qp],3)*_polar_y[_qp] + 8*_polar_x[_qp]*std::pow(_polar_y[_qp],3)) + 
   4*_alpha123*_polar_x[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12*_polar_x[_qp]*_polar_z[_qp] + 4*_alpha123*_polar_x[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp] + 
   _alpha112*(8*std::pow(_polar_x[_qp],3)*_polar_z[_qp] + 8*_polar_x[_qp]*std::pow(_polar_z[_qp],3)));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12*_polar_x[_qp]*_polar_y[_qp] + _alpha112*
    (8*std::pow(_polar_x[_qp],3)*_polar_y[_qp] + 8*_polar_x[_qp]*std::pow(_polar_y[_qp],3)) + 
   4*_alpha123*_polar_x[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12*_polar_y[_qp]*_polar_z[_qp] + 4*_alpha123*std::pow(_polar_x[_qp],2)*_polar_y[_qp]*_polar_z[_qp] + 
   _alpha112*(8*std::pow(_polar_y[_qp],3)*_polar_z[_qp] + 8*_polar_y[_qp]*std::pow(_polar_z[_qp],3)));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12*_polar_x[_qp]*_polar_z[_qp] + 4*_alpha123*_polar_x[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp] + 
   _alpha112*(8*std::pow(_polar_x[_qp],3)*_polar_z[_qp] + 8*_polar_x[_qp]*std::pow(_polar_z[_qp],3)));
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12*_polar_y[_qp]*_polar_z[_qp] + 4*_alpha123*std::pow(_polar_x[_qp],2)*_polar_y[_qp]*_polar_z[_qp] + 
   _alpha112*(8*std::pow(_polar_y[_qp],3)*_polar_z[_qp] + 8*_polar_y[_qp]*std::pow(_polar_z[_qp],3)));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
