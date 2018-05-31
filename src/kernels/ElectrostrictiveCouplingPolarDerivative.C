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

#include "ElectrostrictiveCouplingPolarDerivative.h"

class ElectrostrictiveCouplingPolarDerivative;

registerMooseObject("FerretApp", ElectrostrictiveCouplingPolarDerivative);

template<>
InputParameters validParams<ElectrostrictiveCouplingPolarDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t polarization of the electrostrictive coupling energy. Note: for BFO only.");
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elastic displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<Real>("q11", "the 11 component of rotostrictive coupling tensor");
  params.addRequiredParam<Real>("q12", "the 12 component of rotostrictive coupling tensor");
  params.addRequiredParam<Real>("q44", "the 44 component of rotostrictive coupling tensor");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

ElectrostrictiveCouplingPolarDerivative::ElectrostrictiveCouplingPolarDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _q11(getParam<Real>("q11")),
   _q12(getParam<Real>("q12")),
   _q44(getParam<Real>("q44")),
   _len_scale(getParam<Real>("len_scale"))
{
  std::cout<<"WARNING: This kernel is only for use with BiFeO3 potentials. Its form has not been tested for other materials (like PbTiO3). Did you mean to use FerroelectricCouplingP?"<<"\n";
}

Real
ElectrostrictiveCouplingPolarDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return -_test[_i][_qp] * (-2*_polar_x[_qp]*_q11*_disp_x_grad[_qp](0) - 2*_q44*((_polar_y[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)))/2. + (_polar_z[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))/2.) - _q12*(2*_polar_x[_qp]*_disp_y_grad[_qp](1) + 2*_polar_x[_qp]*_disp_z_grad[_qp](2)));
  }
  else if (_component == 1)
  {
    return -_test[_i][_qp] * (-2*_polar_y[_qp]*_q11*_disp_y_grad[_qp](1) - 2*_q44*((_polar_x[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)))/2. + (_polar_z[_qp]*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))/2.) - _q12*(2*_polar_y[_qp]*_disp_x_grad[_qp](0) + 2*_polar_y[_qp]*_disp_z_grad[_qp](2)));
  }
  else if (_component == 2)
  {
    return -_test[_i][_qp] * (-(_q12*(2*_polar_z[_qp]*_disp_x_grad[_qp](0) + 2*_polar_z[_qp]*_disp_y_grad[_qp](1))) - 2*_q44*((_polar_x[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))/2. + (_polar_y[_qp]*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))/2.) - 2*_polar_z[_qp]*_q11*_disp_z_grad[_qp](2));
  }
  else
    return 0.0;
}

Real
ElectrostrictiveCouplingPolarDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (-2*_q11*_disp_x_grad[_qp](0) - _q12*(2*_disp_y_grad[_qp](1) + 2*_disp_z_grad[_qp](2)));
  }
  else if (_component == 1)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (-2*_q11*_disp_y_grad[_qp](1) - _q12*(2*_disp_x_grad[_qp](0) + 2*_disp_z_grad[_qp](2)));
  }
  else if (_component == 2)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (-(_q12*(2*_disp_x_grad[_qp](0) + 2*_disp_y_grad[_qp](1))) - 2*_q11*_disp_z_grad[_qp](2));
  }
  else
    return 0.0;
}

Real
ElectrostrictiveCouplingPolarDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return -_test[_i][_qp] * _phi[_j][_qp] * (-(_q44*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0))));
    }
    else if (jvar == _polar_z_var)
    {
      return -_test[_i][_qp] * _phi[_j][_qp] * (-(_q44*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0))));
    }
    else if (jvar == _disp_x_var)
    {
      return -_test[_i][_qp] * (-2*_polar_x[_qp]*_grad_phi[_j][_qp](0)*_q11 - 2*((_polar_y[_qp]*_grad_phi[_j][_qp](1))/2. + (_polar_z[_qp]*_grad_phi[_j][_qp](2))/2.)*_q44);
    }
    else if (jvar == _disp_y_var)
    {
      return -_test[_i][_qp] * (-2*_polar_x[_qp]*_grad_phi[_j][_qp](1)*_q12 - _polar_y[_qp]*_grad_phi[_j][_qp](0)*_q44);
    }
    else if (jvar == _disp_z_var)
    {
      return -_test[_i][_qp] * (-2*_polar_x[_qp]*_grad_phi[_j][_qp](2)*_q12 - _polar_z[_qp]*_grad_phi[_j][_qp](0)*_q44);
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
      return -_test[_i][_qp] * (-(_q44*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0))));
    }
    else if (jvar == _polar_z_var)
    {
      return -_test[_i][_qp] * (-(_q44*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1))));
    }
    else if (jvar == _disp_x_var)
    {
      return -_test[_i][_qp] * (-2*_polar_y[_qp]*_grad_phi[_j][_qp](0)*_q12 - _polar_x[_qp]*_grad_phi[_j][_qp](1)*_q44);
    }
    else if (jvar == _disp_y_var)
    {
      return -_test[_i][_qp] * (-2*_polar_y[_qp]*_grad_phi[_j][_qp](1)*_q11 - 2*((_polar_x[_qp]*_grad_phi[_j][_qp](0))/2. + (_polar_z[_qp]*_grad_phi[_j][_qp](2))/2.)*_q44);
    }
    else if (jvar == _disp_z_var)
    {
      return -_test[_i][_qp] * (-2*_polar_y[_qp]*_grad_phi[_j][_qp](2)*_q12 - _polar_z[_qp]*_grad_phi[_j][_qp](1)*_q44);
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
      return -_test[_i][_qp] * (-(_q44*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0))));
    }
    else if (jvar == _polar_y_var)
    {
      return -_test[_i][_qp] * (-(_q44*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0))));
    }
    else if (jvar == _disp_x_var)
    {
      return -_test[_i][_qp] * (-2*_polar_z[_qp]*_grad_phi[_j][_qp](0)*_q12 - _polar_x[_qp]*_grad_phi[_j][_qp](2)*_q44);
    }
    else if (jvar == _disp_y_var)
    {
      return -_test[_i][_qp] * (-2*_polar_z[_qp]*_grad_phi[_j][_qp](1)*_q12 - _polar_y[_qp]*_grad_phi[_j][_qp](2)*_q44);
    }
    else if (jvar == _disp_z_var)
    {
      return -_test[_i][_qp] * (-2*_polar_z[_qp]*_grad_phi[_j][_qp](2)*_q11 - 2*((_polar_x[_qp]*_grad_phi[_j][_qp](0))/2. + (_polar_y[_qp]*_grad_phi[_j][_qp](1))/2.)*_q44);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
