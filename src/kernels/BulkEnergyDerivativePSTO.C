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

#include "BulkEnergyDerivativePSTO.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", BulkEnergyDerivativePSTO);

template<>
InputParameters validParams<BulkEnergyDerivativePSTO>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("a1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("T", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("Tc", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("a2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("a3", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("a4", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("a5", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("eps", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x3", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x4", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x5", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x6", "The coefficients of the Landau expansion");
  return params;
}

BulkEnergyDerivativePSTO::BulkEnergyDerivativePSTO(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _a1(getParam<Real>("a1")),
   _T(getParam<Real>("T")),
   _Tc(getParam<Real>("Tc")),
   _a2(getParam<Real>("a2")),
   _a3(getParam<Real>("a3")),
   _a4(getParam<Real>("a4")),
   _a5(getParam<Real>("a5")),
   _eps(getParam<Real>("eps")),
   _x1(getParam<Real>("x1")),
   _x2(getParam<Real>("x2")),
   _x3(getParam<Real>("x3")),
   _x4(getParam<Real>("x4")),
   _x5(getParam<Real>("x5")),
   _x6(getParam<Real>("x6"))
{
}

Real
BulkEnergyDerivativePSTO::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp]*(4*_a2*Utility::pow<3>(_polar_x[_qp]) + 6*_a4*Utility::pow<5>(_polar_x[_qp]) + 2*_a3*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + _a5*(4*Utility::pow<3>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*_polar_x[_qp]*Utility::pow<4>(_polar_y[_qp])) + 2*_a1*_polar_x[_qp]*(_T - _Tc) + _eps*(2*_polar_x[_qp]*_x1 + 4*Utility::pow<3>(_polar_x[_qp])*_x2 + 2*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_x3) + 
   Utility::pow<2>(_eps)*(2*_polar_x[_qp]*_x4 + 4*Utility::pow<3>(_polar_x[_qp])*_x5 + 2*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_x6));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*(2*_a3*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 4*_a2*Utility::pow<3>(_polar_y[_qp]) + 6*_a4*Utility::pow<5>(_polar_y[_qp]) + _a5*(2*Utility::pow<4>(_polar_x[_qp])*_polar_y[_qp] + 4*Utility::pow<2>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp])) + 2*_a1*_polar_y[_qp]*(_T - _Tc) + _eps*(2*_polar_y[_qp]*_x1 + 4*Utility::pow<3>(_polar_y[_qp])*_x2 + 2*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_x3) + 
   Utility::pow<2>(_eps)*(2*_polar_y[_qp]*_x4 + 4*Utility::pow<3>(_polar_y[_qp])*_x5 + 2*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_x6));
  }
  else if (_component == 2)
  {
    return 0.0;
  }
  else
    return 0.0;
}

Real
BulkEnergyDerivativePSTO::computeQpJacobian()
{
  if (_component == 0)
  {
    return _phi[_j][_qp]*_test[_i][_qp]*(12*_a2*Utility::pow<2>(_polar_x[_qp]) + 30*_a4*Utility::pow<4>(_polar_x[_qp]) + 2*_a3*Utility::pow<2>(_polar_y[_qp]) + _a5*(12*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<4>(_polar_y[_qp])) + 2*_a1*(_T - _Tc) + _eps*(2*_x1 + 12*Utility::pow<2>(_polar_x[_qp])*_x2 + 2*Utility::pow<2>(_polar_y[_qp])*_x3) + 
   Utility::pow<2>(_eps)*(2*_x4 + 12*Utility::pow<2>(_polar_x[_qp])*_x5 + 2*Utility::pow<2>(_polar_y[_qp])*_x6));
  }
  else if (_component == 1)
  {
    return _phi[_j][_qp]*_test[_i][_qp]*(2*_a3*Utility::pow<2>(_polar_x[_qp]) + 12*_a2*Utility::pow<2>(_polar_y[_qp]) + 30*_a4*Utility::pow<4>(_polar_y[_qp]) + _a5*(2*Utility::pow<4>(_polar_x[_qp]) + 12*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])) + 2*_a1*(_T - _Tc) + _eps*(2*_x1 + 12*Utility::pow<2>(_polar_y[_qp])*_x2 + 2*Utility::pow<2>(_polar_x[_qp])*_x3) + 
   Utility::pow<2>(_eps)*(2*_x4 + 12*Utility::pow<2>(_polar_y[_qp])*_x5 + 2*Utility::pow<2>(_polar_x[_qp])*_x6));
  }
  else if (_component == 2)
  {
    return 0.0;
  }
  else
    return 0.0;
}

Real
BulkEnergyDerivativePSTO::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp]*_test[_i][_qp]*(4*_a3*_polar_x[_qp]*_polar_y[_qp] + _a5*(8*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + 8*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp])) + 4*_eps*_polar_x[_qp]*_polar_y[_qp]*_x3 + 4*Utility::pow<2>(_eps)*_polar_x[_qp]*_polar_y[_qp]*_x6);
    }
    else if (jvar == _polar_z_var)
    {
      return 0.0;
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
      return _phi[_j][_qp]*_test[_i][_qp]*(4*_a3*_polar_x[_qp]*_polar_y[_qp] + _a5*(8*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + 8*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp])) + 4*_eps*_polar_x[_qp]*_polar_y[_qp]*_x3 + 4*Utility::pow<2>(_eps)*_polar_x[_qp]*_polar_y[_qp]*_x6);
    }
    else if (jvar == _polar_z_var)
    {
      return 0.0;
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
      return 0.0;
    }
    else if (jvar == _polar_y_var)
    {
      return 0.0;
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
