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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "BulkEnergyDerivativeSixthCoupledT.h"
#include<cmath>

registerMooseObject("FerretApp", BulkEnergyDerivativeSixthCoupledT);

InputParameters BulkEnergyDerivativeSixthCoupledT::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual for the local free energy which is an sixth order expansion in the "
                             "polarization coupled to the thermal field through the first Landau coefficient.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("temperature", "The temperature at the grid point");
  params.addRequiredParam<Real>("alpha0", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("Tc", "Transition temperature of unstrained ferroelectric");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

BulkEnergyDerivativeSixthCoupledT::BulkEnergyDerivativeSixthCoupledT(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _temperature_var(coupled("temperature")),
   _temperature(coupledValue("temperature")),
   _alpha0(getParam<Real>("alpha0")),
   _alpha11(getParam<Real>("alpha11")),
   _alpha12(getParam<Real>("alpha12")),
   _alpha111(getParam<Real>("alpha111")),
   _alpha112(getParam<Real>("alpha112")),
   _alpha123(getParam<Real>("alpha123")),
   _Tc(getParam<Real>("Tc")),
   _len_scale(getParam<Real>("len_scale"))
{
  std::cout<<"_alpha0 = "<<_alpha0<<"\n";
  std::cout<<"_alpha11 = "<<_alpha11<<"\n";
  std::cout<<"_alpha12 = "<<_alpha12<<"\n";
  std::cout<<"_alpha111 = "<<_alpha111<<"\n";
  std::cout<<"_alpha112 = "<<_alpha112<<"\n";
  std::cout<<"_alpha123 = "<<_alpha123<<"\n";
}

Real
BulkEnergyDerivativeSixthCoupledT::computeQpResidual()
{
  const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1) ? _polar_y: _polar_z;
  const VariableValue & _polar_j = (_component == 0) ? _polar_y : (_component == 1) ? _polar_z: _polar_x;
  const VariableValue & _polar_k = (_component == 0) ? _polar_z : (_component == 1) ? _polar_x: _polar_y;
  Real Rbulk = 0.0;
  Rbulk += ((2.0 * _alpha0 * (_temperature[_qp] - _Tc) * _polar_i[_qp] + 4 * _alpha11 * Utility::pow<3>(_polar_i[_qp]) + 2.0 * _alpha12 * _polar_i[_qp]*(Utility::pow<2>(_polar_j[_qp]) + Utility::pow<2>(_polar_k[_qp])) +
	  6.0 * _alpha111 * Utility::pow<5>(_polar_i[_qp]) + 4.0 * _alpha112 * Utility::pow<3>(_polar_i[_qp]) * (_polar_j[_qp] * _polar_j[_qp]+_polar_k[_qp] * _polar_k[_qp]) +
	  2.0 * _alpha112 * _polar_i[_qp]*(Utility::pow<4>(_polar_j[_qp]) + Utility::pow<4>(_polar_k[_qp])) + 2.0 * _alpha123 * _polar_i[_qp]*Utility::pow<2>(_polar_j[_qp]) * Utility::pow<2>(_polar_k[_qp])) * _test[_i][_qp]) * std::pow(_len_scale, 3.0);
  ///  Moose::out << "\n R_bulk-"; std::cout << _component << " = " << Rbulk;
  return Rbulk;
}

Real
BulkEnergyDerivativeSixthCoupledT::computeQpJacobian()
{
  const VariableValue & _polar_i = (_component == 0)? _polar_x : (_component == 1)? _polar_y: _polar_z;
  const VariableValue & _polar_j = (_component == 0)? _polar_y : (_component == 1)? _polar_z: _polar_x;
  const VariableValue & _polar_k = (_component == 0)? _polar_z : (_component == 1)? _polar_x: _polar_y;
  return (2 * _alpha0 * (_temperature[_qp] - _Tc) + 12.0 * _alpha11 * Utility::pow<2>(_polar_i[_qp]) +
	  2.0 * _alpha12 * (Utility::pow<2>(_polar_j[_qp]) + Utility::pow<2>(_polar_k[_qp])) + 30.0 * _alpha111 * Utility::pow<4>(_polar_i[_qp]) +
	  12.0 * _alpha112 * Utility::pow<2>(_polar_i[_qp]) * (Utility::pow<2>(_polar_j[_qp]) + Utility::pow<2>(_polar_k[_qp])) + 2.0 * _alpha112 * (Utility::pow<4>(_polar_j[_qp]) + Utility::pow<4>(_polar_k[_qp])) +
	  2.0 * _alpha123 * Utility::pow<2>(_polar_j[_qp]) * Utility::pow<2>(_polar_k[_qp])
  ) * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
}

Real
BulkEnergyDerivativeSixthCoupledT::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real r;
  //mooseAssert(jvar != variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var)
    {
      const VariableValue & _polar_i = (_component == 0)? _polar_x : (_component == 1)? _polar_y: _polar_z;
      const VariableValue & _polar_j = (jvar == _polar_x_var)? _polar_x : (jvar == _polar_y_var)? _polar_y: _polar_z;
      const VariableValue & _polar_k = ((_component == 0 && jvar == _polar_y_var) || (_component == 1 && jvar == _polar_x_var) )? _polar_z : ( (_component == 0 && jvar == _polar_z_var) || (_component == 2 && jvar == _polar_x_var))? _polar_y: _polar_x;
      r = (4.0 * _alpha12 * _polar_i[_qp] * _polar_j[_qp] + 8.0 * _alpha112 * Utility::pow<3>(_polar_i[_qp]) * _polar_j[_qp]
      + 8.0 *_alpha112 * _polar_i[_qp] * Utility::pow<3>(_polar_j[_qp]) + 4.0 * _alpha123 * _polar_i[_qp] * _polar_j[_qp] * Utility::pow<2>(_polar_k[_qp]));
      return r * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
    }
  if(jvar == _temperature_var)
  {
    const VariableValue& _polar_i = (_component == 0)? _polar_x : (_component == 1)? _polar_y: _polar_z;
    r = 2.0 * _alpha0 * _polar_i[_qp];
    return r * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
  }
  else
    return 0.0;
}
