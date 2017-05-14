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

#include "BulkEnergyDerivativeSixthCoupledT.h"
#include<cmath>

template<>
InputParameters validParams<BulkEnergyDerivativeSixthCoupledT>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
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
  Rbulk += ((2.0 * _alpha0 * (_temperature[_qp] - _Tc) * _polar_i[_qp] + 4 * _alpha11 * std::pow(_polar_i[_qp], 3.0) + 2.0 * _alpha12 * _polar_i[_qp]*(std::pow(_polar_j[_qp], 2.0) + std::pow(_polar_k[_qp], 2.0)) +
	  6.0 * _alpha111 * std::pow(_polar_i[_qp], 5.0) + 4.0 * _alpha112 * std::pow(_polar_i[_qp], 3.0) * (_polar_j[_qp] * _polar_j[_qp]+_polar_k[_qp] * _polar_k[_qp]) +
	  2.0 * _alpha112 * _polar_i[_qp]*(std::pow(_polar_j[_qp], 4.0) + std::pow(_polar_k[_qp], 4.0)) + 2.0 * _alpha123 * _polar_i[_qp]*std::pow(_polar_j[_qp], 2.0) * std::pow(_polar_k[_qp], 2.0)) * _test[_i][_qp]) * std::pow(_len_scale, 3.0);
  ///  Moose::out << "\n R_bulk-"; std::cout << _component << " = " << Rbulk;
  return Rbulk;
}

Real
BulkEnergyDerivativeSixthCoupledT::computeQpJacobian()
{
  const VariableValue & _polar_i = (_component == 0)? _polar_x : (_component == 1)? _polar_y: _polar_z;
  const VariableValue & _polar_j = (_component == 0)? _polar_y : (_component == 1)? _polar_z: _polar_x;
  const VariableValue & _polar_k = (_component == 0)? _polar_z : (_component == 1)? _polar_x: _polar_y;
  return (2 * _alpha0 * (_temperature[_qp] - _Tc) + 12.0 * _alpha11 * std::pow(_polar_i[_qp], 2.0) +
	  2.0 * _alpha12 * (std::pow(_polar_j[_qp], 2.0) + std::pow(_polar_k[_qp], 2.0)) + 30.0 * _alpha111 * std::pow(_polar_i[_qp], 4.0) +
	  12.0 * _alpha112 * std::pow(_polar_i[_qp], 2.0) * (std::pow(_polar_j[_qp],2.0) + std::pow(_polar_k[_qp], 2.0)) + 2.0 * _alpha112 * (std::pow(_polar_j[_qp], 4.0) + std::pow(_polar_k[_qp], 4.0)) +
	  2.0 * _alpha123 * std::pow(_polar_j[_qp], 2.0) * std::pow(_polar_k[_qp], 2.0)
  ) * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
}

Real
BulkEnergyDerivativeSixthCoupledT::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real r;
  mooseAssert(jvar != variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var)
    {
      const VariableValue & _polar_i = (_component == 0)? _polar_x : (_component == 1)? _polar_y: _polar_z;
      const VariableValue & _polar_j = (jvar == _polar_x_var)? _polar_x : (jvar == _polar_y_var)? _polar_y: _polar_z;
      const VariableValue & _polar_k = ((_component == 0 && jvar == _polar_y_var) || (_component == 1 && jvar == _polar_x_var) )? _polar_z : ( (_component == 0 && jvar == _polar_z_var) || (_component == 2 && jvar == _polar_x_var))? _polar_y: _polar_x;
      r = (4.0 * _alpha12 * _polar_i[_qp] * _polar_j[_qp] + 8.0 * _alpha112 * std::pow(_polar_i[_qp], 3.0) * _polar_j[_qp]
      + 8.0 *_alpha112 * _polar_i[_qp] * std::pow(_polar_j[_qp], 3.0) + 4.0 * _alpha123 * _polar_i[_qp] * _polar_j[_qp] * std::pow(_polar_k[_qp], 2.0));
      return r * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
    }
  if(jvar == _temperature_var)
  {
    const VariableValue& _polar_i= (_component == 0)? _polar_x : (_component == 1)? _polar_y: _polar_z;
    const VariableValue& _polar_j= (_component == 0)? _polar_y : (_component == 1)? _polar_z: _polar_x;
    const VariableValue& _polar_k= (_component == 0)? _polar_z : (_component == 1)? _polar_x: _polar_y;
    r = 2.0 * _alpha0 * _polar_i[_qp];
    return r * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
  }
  else
    return 0.0;
}
