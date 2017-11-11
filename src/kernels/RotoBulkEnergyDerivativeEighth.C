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

#include "RotoBulkEnergyDerivativeEighth.h"
#include<cmath>

template<>
InputParameters validParams<RotoBulkEnergyDerivativeEighth>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_z", "The z component of the antiferrodistortive tilt");
  params.addRequiredParam<Real>("beta1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta123", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1122", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotoBulkEnergyDerivativeEighth::RotoBulkEnergyDerivativeEighth(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _beta1(getParam<Real>("beta1")),
   _beta11(getParam<Real>("beta11")),
   _beta12(getParam<Real>("beta12")),
   _beta111(getParam<Real>("beta111")),
   _beta112(getParam<Real>("beta112")),
   _beta123(getParam<Real>("beta123")),
   _beta1111(getParam<Real>("beta1111")),
   _beta1112(getParam<Real>("beta1112")),
   _beta1122(getParam<Real>("beta1122")),
   _beta1123(getParam<Real>("beta1123")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotoBulkEnergyDerivativeEighth::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_beta1*_antiferrodis_A_x[_qp] + 4*_beta11*std::pow(_antiferrodis_A_x[_qp],3) + 6*_beta111*std::pow(_antiferrodis_A_x[_qp],5) + 8*_beta1111*std::pow(_antiferrodis_A_x[_qp],7) + (_beta123*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2))/3. + (_beta12*(2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2) + 2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)))/2. + 
   (_beta1122*(4*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_antiferrodis_A_y[_qp],4) + 4*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_antiferrodis_A_z[_qp],4)))/6. + (_beta1123*(6*std::pow(_antiferrodis_A_x[_qp],5)*std::pow(_antiferrodis_A_z[_qp],2) + 2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],4)*std::pow(_antiferrodis_A_z[_qp],2) + 2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],4)))/12. + 
   (_beta112*(2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],4) + 2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],4) + 4*std::pow(_antiferrodis_A_x[_qp],3)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/3. + (_beta1112*(2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],6) + 2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],6) + 6*std::pow(_antiferrodis_A_x[_qp],5)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/4.);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*_beta1*_antiferrodis_A_y[_qp] + 4*_beta11*std::pow(_antiferrodis_A_y[_qp],3) + 6*_beta111*std::pow(_antiferrodis_A_y[_qp],5) + 8*_beta1111*std::pow(_antiferrodis_A_y[_qp],7) + (_beta123*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2))/3. + (_beta12*(2*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)))/2. + 
   (_beta1123*(4*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_antiferrodis_A_z[_qp],2) + 2*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],4)))/12. + (_beta1122*(4*std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_antiferrodis_A_y[_qp],3) + 4*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_antiferrodis_A_z[_qp],4)))/6. + 
   (_beta112*(2*std::pow(_antiferrodis_A_x[_qp],4)*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],4) + 4*std::pow(_antiferrodis_A_y[_qp],3)*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/3. + (_beta1112*(2*std::pow(_antiferrodis_A_x[_qp],6)*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],6) + 6*std::pow(_antiferrodis_A_y[_qp],5)*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/4.);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_beta1*_antiferrodis_A_z[_qp] + (_beta123*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp])/3. + 4*_beta11*std::pow(_antiferrodis_A_z[_qp],3) + 6*_beta111*std::pow(_antiferrodis_A_z[_qp],5) + 8*_beta1111*std::pow(_antiferrodis_A_z[_qp],7) + (_beta12*(2*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp] + 2*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]))/2. + 
   (_beta1123*(2*std::pow(_antiferrodis_A_x[_qp],6)*_antiferrodis_A_z[_qp] + 2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],4)*_antiferrodis_A_z[_qp] + 4*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],3)))/12. + (_beta1122*(4*std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_antiferrodis_A_z[_qp],3) + 4*std::pow(_antiferrodis_A_y[_qp],4)*std::pow(_antiferrodis_A_z[_qp],3)))/6. + 
   (_beta112*(2*std::pow(_antiferrodis_A_x[_qp],4)*_antiferrodis_A_z[_qp] + 2*std::pow(_antiferrodis_A_y[_qp],4)*_antiferrodis_A_z[_qp] + 4*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],3)))/3. + (_beta1112*(2*std::pow(_antiferrodis_A_x[_qp],6)*_antiferrodis_A_z[_qp] + 2*std::pow(_antiferrodis_A_y[_qp],6)*_antiferrodis_A_z[_qp] + 6*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],5)))/4.);
  }
  else
    return 0.0;
}

Real
RotoBulkEnergyDerivativeEighth::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1 + 12*_beta11*std::pow(_antiferrodis_A_x[_qp],2) + 30*_beta111*std::pow(_antiferrodis_A_x[_qp],4) + 56*_beta1111*std::pow(_antiferrodis_A_x[_qp],6) + (_beta123*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2))/3. + (_beta12*(2*std::pow(_antiferrodis_A_y[_qp],2) + 2*std::pow(_antiferrodis_A_z[_qp],2)))/2. + 
   (_beta1122*(12*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],4) + 12*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_z[_qp],4)))/6. + (_beta1123*(30*std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_antiferrodis_A_z[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],4)*std::pow(_antiferrodis_A_z[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],4)))/12. + 
   (_beta112*(2*std::pow(_antiferrodis_A_y[_qp],4) + 2*std::pow(_antiferrodis_A_z[_qp],4) + 12*std::pow(_antiferrodis_A_x[_qp],2)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/3. + (_beta1112*(2*std::pow(_antiferrodis_A_y[_qp],6) + 2*std::pow(_antiferrodis_A_z[_qp],6) + 30*std::pow(_antiferrodis_A_x[_qp],4)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/4.);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1 + 12*_beta11*std::pow(_antiferrodis_A_y[_qp],2) + 30*_beta111*std::pow(_antiferrodis_A_y[_qp],4) + 56*_beta1111*std::pow(_antiferrodis_A_y[_qp],6) + (_beta123*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2))/3. + (_beta12*(2*std::pow(_antiferrodis_A_x[_qp],2) + 2*std::pow(_antiferrodis_A_z[_qp],2)))/2. + 
   (_beta1123*(12*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2) + 2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_z[_qp],4)))/12. + (_beta1122*(12*std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_antiferrodis_A_y[_qp],2) + 12*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],4)))/6. + 
   (_beta112*(2*std::pow(_antiferrodis_A_x[_qp],4) + 2*std::pow(_antiferrodis_A_z[_qp],4) + 12*std::pow(_antiferrodis_A_y[_qp],2)*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/3. + (_beta1112*(2*std::pow(_antiferrodis_A_x[_qp],6) + 2*std::pow(_antiferrodis_A_z[_qp],6) + 30*std::pow(_antiferrodis_A_y[_qp],4)*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))))/4.);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1 + (_beta123*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],2))/3. + (_beta12*(2*std::pow(_antiferrodis_A_x[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],2)))/2. + 12*_beta11*std::pow(_antiferrodis_A_z[_qp],2) + 30*_beta111*std::pow(_antiferrodis_A_z[_qp],4) + 56*_beta1111*std::pow(_antiferrodis_A_z[_qp],6) + 
   (_beta1123*(2*std::pow(_antiferrodis_A_x[_qp],6) + 2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],4) + 12*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2)))/12. + (_beta1122*(12*std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_antiferrodis_A_z[_qp],2) + 12*std::pow(_antiferrodis_A_y[_qp],4)*std::pow(_antiferrodis_A_z[_qp],2)))/6. + 
   (_beta112*(2*std::pow(_antiferrodis_A_x[_qp],4) + 2*std::pow(_antiferrodis_A_y[_qp],4) + 12*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],2)))/3. + (_beta1112*(2*std::pow(_antiferrodis_A_x[_qp],6) + 2*std::pow(_antiferrodis_A_y[_qp],6) + 30*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],4)))/4.);
  }
  else
    return 0.0;
}

Real
RotoBulkEnergyDerivativeEighth::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (2*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp] + (8*_beta1122*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_antiferrodis_A_y[_qp],3))/3. + (_beta112*(8*std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_y[_qp] + 8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],3)))/3. + (_beta1112*(12*std::pow(_antiferrodis_A_x[_qp],5)*_antiferrodis_A_y[_qp] + 12*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],5)))/4. + (2*_beta123*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2))/3. + (_beta1123*(8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_antiferrodis_A_z[_qp],2) + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],4)))/12.);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (2*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp] + (2*_beta123*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp])/3. + (8*_beta1122*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_antiferrodis_A_z[_qp],3))/3. + (_beta112*(8*std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],3)))/3. + 
   (_beta1123*(12*std::pow(_antiferrodis_A_x[_qp],5)*_antiferrodis_A_z[_qp] + 4*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],4)*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],3)))/12. + (_beta1112*(12*std::pow(_antiferrodis_A_x[_qp],5)*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],5)))/4.);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (2*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp] + (8*_beta1122*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_antiferrodis_A_y[_qp],3))/3. + (_beta112*(8*std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_y[_qp] + 8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],3)))/3. + (_beta1112*(12*std::pow(_antiferrodis_A_x[_qp],5)*_antiferrodis_A_y[_qp] + 12*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],5)))/4. + (2*_beta123*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2))/3. + (_beta1123*(8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_antiferrodis_A_z[_qp],2) + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],4)))/12.);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (2*_beta12*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + (2*_beta123*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp])/3. + (8*_beta1122*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_antiferrodis_A_z[_qp],3))/3. + (_beta112*(8*std::pow(_antiferrodis_A_y[_qp],3)*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],3)))/3. + (_beta1123*(8*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],3)*_antiferrodis_A_z[_qp] + 8*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],3)))/12. + (_beta1112*(12*std::pow(_antiferrodis_A_y[_qp],5)*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],5)))/4.);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (2*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp] + (2*_beta123*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp])/3. + (8*_beta1122*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_antiferrodis_A_z[_qp],3))/3. + (_beta112*(8*std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],3)))/3. + 
   (_beta1123*(12*std::pow(_antiferrodis_A_x[_qp],5)*_antiferrodis_A_z[_qp] + 4*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],4)*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],3)))/12. + (_beta1112*(12*std::pow(_antiferrodis_A_x[_qp],5)*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],5)))/4.);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (2*_beta12*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + (2*_beta123*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp])/3. + (8*_beta1122*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_antiferrodis_A_z[_qp],3))/3. + (_beta112*(8*std::pow(_antiferrodis_A_y[_qp],3)*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],3)))/3. + (_beta1123*(8*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],3)*_antiferrodis_A_z[_qp] + 8*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],3)))/12. + (_beta1112*(12*std::pow(_antiferrodis_A_y[_qp],5)*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],5)))/4.);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
