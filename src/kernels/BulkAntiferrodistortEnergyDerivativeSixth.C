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

#include "BulkAntiferrodistortEnergyDerivativeSixth.h"
#include "libmesh/utility.h"
#include<cmath>

template<>
InputParameters validParams<BulkAntiferrodistortEnergyDerivativeSixth>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_y", 0.0, "The y component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt vector");
  params.addRequiredParam<Real>("beta1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

BulkAntiferrodistortEnergyDerivativeSixth::BulkAntiferrodistortEnergyDerivativeSixth(const InputParameters & parameters)
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
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkAntiferrodistortEnergyDerivativeSixth::computeQpResidual()
{
  const VariableValue & _antiferrodis_A_i = (_component == 0) ? _antiferrodis_A_x : (_component == 1) ? _antiferrodis_A_y: _antiferrodis_A_z;
  const VariableValue & _antiferrodis_A_j = (_component == 0) ? _antiferrodis_A_y : (_component == 1) ? _antiferrodis_A_z: _antiferrodis_A_x;
  const VariableValue & _antiferrodis_A_k = (_component == 0) ? _antiferrodis_A_z : (_component == 1) ? _antiferrodis_A_x: _antiferrodis_A_y;
  Real Rbulk = 0.0;
  Rbulk += ((2.0 * _beta1 * _antiferrodis_A_i[_qp] + 4.0 * _beta11 * Utility::pow<3>(_antiferrodis_A_i[_qp]) + 2.0 * _beta12 * _antiferrodis_A_i[_qp]*(Utility::pow<2>(_antiferrodis_A_j[_qp]) + Utility::pow<2>(_antiferrodis_A_k[_qp])) +
	  6.0 * _beta111 * std::pow(_antiferrodis_A_i[_qp], 5.0) + 4.0 * _beta112 * Utility::pow<3>(_antiferrodis_A_i[_qp]) * (_antiferrodis_A_j[_qp] * _antiferrodis_A_j[_qp]+_antiferrodis_A_k[_qp] * _antiferrodis_A_k[_qp]) +
	  2.0 * _beta112 * _antiferrodis_A_i[_qp]*(Utility::pow<4>(_antiferrodis_A_j[_qp]) + Utility::pow<4>(_antiferrodis_A_k[_qp])) + 2.0 * _beta123 * _antiferrodis_A_i[_qp]*Utility::pow<2>(_antiferrodis_A_j[_qp]) * Utility::pow<2>(_antiferrodis_A_k[_qp])) * _test[_i][_qp]) * Utility::pow<3>(_len_scale);
  ///  Moose::out << "\n R_bulk-"; std::cout << _component << " = " << Rbulk;
  return Rbulk;
}

Real
BulkAntiferrodistortEnergyDerivativeSixth::computeQpJacobian()
{
  const VariableValue & _antiferrodis_A_i = (_component == 0)? _antiferrodis_A_x : (_component == 1)? _antiferrodis_A_y: _antiferrodis_A_z;
  const VariableValue & _antiferrodis_A_j = (_component == 0)? _antiferrodis_A_y : (_component == 1)? _antiferrodis_A_z: _antiferrodis_A_x;
  const VariableValue & _antiferrodis_A_k = (_component == 0)? _antiferrodis_A_z : (_component == 1)? _antiferrodis_A_x: _antiferrodis_A_y;
  return (2.0 * _beta1 + 12.0 * _beta11 * std::pow(_antiferrodis_A_i[_qp], 2) +
	  2.0 * _beta12 * (Utility::pow<2>(_antiferrodis_A_j[_qp]) + Utility::pow<2>(_antiferrodis_A_k[_qp])) + 30.0 * _beta111 * Utility::pow<4>(_antiferrodis_A_i[_qp]) +
	  12.0 * _beta112 * Utility::pow<2>(_antiferrodis_A_i[_qp]) * (Utility::pow<2>(_antiferrodis_A_j[_qp]) + Utility::pow<2>(_antiferrodis_A_k[_qp])) + 2.0 * _beta112 * (Utility::pow<4>(_antiferrodis_A_j[_qp]) + Utility::pow<4>(_antiferrodis_A_k[_qp])) +
	  2.0 * _beta123 * Utility::pow<2>(_antiferrodis_A_j[_qp]) * Utility::pow<2>(_antiferrodis_A_k[_qp])
  ) * _test[_i][_qp] * _phi[_j][_qp] * Utility::pow<3>(_len_scale);
}

Real
BulkAntiferrodistortEnergyDerivativeSixth::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real r;
  //mooseAssert(jvar != variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_antiferrodis_A_x_var || jvar==_antiferrodis_A_y_var || jvar==_antiferrodis_A_z_var)
    {
      const VariableValue & _antiferrodis_A_i = (_component == 0)? _antiferrodis_A_x : (_component == 1)? _antiferrodis_A_y: _antiferrodis_A_z;
      const VariableValue & _antiferrodis_A_j = (jvar == _antiferrodis_A_x_var)? _antiferrodis_A_x : (jvar == _antiferrodis_A_y_var)? _antiferrodis_A_y: _antiferrodis_A_z;
      const VariableValue & _antiferrodis_A_k = ((_component == 0 && jvar == _antiferrodis_A_y_var) || (_component == 1 && jvar == _antiferrodis_A_x_var) )? _antiferrodis_A_z : ( (_component == 0 && jvar == _antiferrodis_A_z_var) || (_component == 2 && jvar == _antiferrodis_A_x_var))? _antiferrodis_A_y: _antiferrodis_A_x;
      r = (4.0 * _beta12 * _antiferrodis_A_i[_qp] * _antiferrodis_A_j[_qp] + 8.0 * _beta112 * Utility::pow<3>(_antiferrodis_A_i[_qp]) * _antiferrodis_A_j[_qp]
      + 8.0 *_beta112 * _antiferrodis_A_i[_qp] * Utility::pow<3>(_antiferrodis_A_j[_qp]) + 4.0 * _beta123 * _antiferrodis_A_i[_qp] * _antiferrodis_A_j[_qp] * Utility::pow<2>(_antiferrodis_A_k[_qp]));
      return r * _test[_i][_qp] * _phi[_j][_qp] * Utility::pow<3>(_len_scale);
    }
  else
    return 0.0;
}
