/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3.0 of the License, or
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

#include "RotoBulkEnergyDerivativeEighthAlt.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotoBulkEnergyDerivativeEighthAlt);

template<>
InputParameters validParams<RotoBulkEnergyDerivativeEighthAlt>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2.0 for z)");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt");
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

RotoBulkEnergyDerivativeEighthAlt::RotoBulkEnergyDerivativeEighthAlt(const InputParameters & parameters)
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
RotoBulkEnergyDerivativeEighthAlt::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_beta1*_antiferrodis_A_x[_qp] + 4*_beta11*Utility::pow<3>(_antiferrodis_A_x[_qp]) + 6*_beta111*Utility::pow<5>(_antiferrodis_A_x[_qp]) + 8*_beta1111*Utility::pow<7>(_antiferrodis_A_x[_qp]) + 
   2*_beta123*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12*(2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1122*(4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112*(2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      4*Utility::pow<3>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112*(2*_antiferrodis_A_x[_qp]*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      6*Utility::pow<5>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*_beta1*_antiferrodis_A_y[_qp] + 4*_beta11*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 6*_beta111*Utility::pow<5>(_antiferrodis_A_y[_qp]) + 8*_beta1111*Utility::pow<7>(_antiferrodis_A_y[_qp]) + 
   2*_beta123*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12*(2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 4*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1122*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      4*Utility::pow<3>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112*(2*Utility::pow<6>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      6*Utility::pow<5>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_beta1*_antiferrodis_A_z[_qp] + 2*_beta123*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 4*_beta11*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   6*_beta111*Utility::pow<5>(_antiferrodis_A_z[_qp]) + 8*_beta1111*Utility::pow<7>(_antiferrodis_A_z[_qp]) + 
   _beta12*(2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]) + 
   _beta1123*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      4*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1122*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 4*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta112*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      4*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1112*(2*Utility::pow<6>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<6>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      6*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<5>(_antiferrodis_A_z[_qp])));
  }
  else
    return 0.0;
}

Real
RotoBulkEnergyDerivativeEighthAlt::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1 + 12*_beta11*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 30*_beta111*Utility::pow<4>(_antiferrodis_A_x[_qp]) + 56*_beta1111*Utility::pow<6>(_antiferrodis_A_x[_qp]) + 
   2*_beta123*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12*(2*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1122*(12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112*(2*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      12*Utility::pow<2>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112*(2*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      30*Utility::pow<4>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1 + 12*_beta11*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 30*_beta111*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 56*_beta1111*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 
   2*_beta123*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12*(2*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1122*(12*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112*(2*Utility::pow<4>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      12*Utility::pow<2>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112*(2*Utility::pow<6>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      30*Utility::pow<4>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1 + 2*_beta123*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp]) + _beta12*(2*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])) + 
   12*_beta11*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 30*_beta111*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 56*_beta1111*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
   _beta1123*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 
      12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1122*(12*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 12*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta112*(2*Utility::pow<4>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 
      12*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1112*(2*Utility::pow<6>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 
      30*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<4>(_antiferrodis_A_z[_qp])));
  }
  else
    return 0.0;
}

Real
RotoBulkEnergyDerivativeEighthAlt::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp] + 16*_beta1122*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 
   _beta112*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])) + 
   _beta1112*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_y[_qp])) + 4*_beta123*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
   _beta1123*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp])));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 16*_beta1122*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 4*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
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
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp] + 16*_beta1122*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 
   _beta112*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])) + 
   _beta1112*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_y[_qp])) + 4*_beta123*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
   _beta1123*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp])));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 16*_beta1122*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112*(8*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 8*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112*(12*Utility::pow<5>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_y[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
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
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 16*_beta1122*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 4*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 16*_beta1122*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112*(8*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 8*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112*(12*Utility::pow<5>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_y[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
