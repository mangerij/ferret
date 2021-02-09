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

#include "InducedPElectrostrictiveCouplingPolarDerivative.h"

class InducedPElectrostrictiveCouplingPolarDerivative;

registerMooseObject("FerretApp", InducedPElectrostrictiveCouplingPolarDerivative);

template<>
InputParameters validParams<InducedPElectrostrictiveCouplingPolarDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t polarization of the electrostrictive coupling energy. Note: for BFO only.");
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elastic displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("induced_polar_x", "The x component of the induced_polarization");
  params.addRequiredCoupledVar("induced_polar_y", "The y component of the induced_polarization");
  params.addCoupledVar("induced_polar_z", 0.0, "The z component of the induced_polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  return params;
}

InducedPElectrostrictiveCouplingPolarDerivative::InducedPElectrostrictiveCouplingPolarDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _induced_polar_x_var(coupled("induced_polar_x")),
   _induced_polar_y_var(coupled("induced_polar_y")),
   _induced_polar_z_var(coupled("induced_polar_z")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _induced_polar_x(coupledValue("induced_polar_x")),
   _induced_polar_y(coupledValue("induced_polar_y")),
   _induced_polar_z(coupledValue("induced_polar_z")),
   _q11(getMaterialProperty<Real>("q11")),
   _q12(getMaterialProperty<Real>("q12")),
   _q44(getMaterialProperty<Real>("q44"))
{
}

Real
InducedPElectrostrictiveCouplingPolarDerivative::computeQpResidual()
{
//as with ElectrostrictiveCouplingPolarDerivative, there is an overall minus sign included here.
  if (_component == 0)
  {
    return -_test[_i][_qp] * ((-(_q44[_qp]*(_induced_polar_y[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)) + _polar_y[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)) + (_induced_polar_z[_qp] + _polar_z[_qp])*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))) - 4*_induced_polar_x[_qp]*(_q11[_qp]*_disp_x_grad[_qp](0) + _q12[_qp]*(_disp_y_grad[_qp](1) + _disp_z_grad[_qp](2))) - 4*_polar_x[_qp]*(_q11[_qp]*_disp_x_grad[_qp](0) + _q12[_qp]*(_disp_y_grad[_qp](1) + _disp_z_grad[_qp](2))))/2.);
  }
  else if (_component == 1)
  {
    return -_test[_i][_qp] * ((-(_q44[_qp]*(_induced_polar_x[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)) + _polar_x[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)) + (_induced_polar_z[_qp] + _polar_z[_qp])*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))) - 4*_induced_polar_y[_qp]*(_q11[_qp]*_disp_y_grad[_qp](1) + _q12[_qp]*(_disp_x_grad[_qp](0) + _disp_z_grad[_qp](2))) - 4*_polar_y[_qp]*(_q11[_qp]*_disp_y_grad[_qp](1) + _q12[_qp]*(_disp_x_grad[_qp](0) + _disp_z_grad[_qp](2))))/2.);
  }
  else if (_component == 2)
  {
    return -_test[_i][_qp] * ((-(_q44[_qp]*(_induced_polar_x[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)) + _polar_x[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)) + (_induced_polar_y[_qp] + _polar_y[_qp])*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))) - 4*_induced_polar_z[_qp]*(_q12[_qp]*(_disp_x_grad[_qp](0) + _disp_y_grad[_qp](1)) + _q11[_qp]*_disp_z_grad[_qp](2)) - 4*_polar_z[_qp]*(_q12[_qp]*(_disp_x_grad[_qp](0) + _disp_y_grad[_qp](1)) + _q11[_qp]*_disp_z_grad[_qp](2)))/2.);
  }
  else
    return 0.0;
}

Real
InducedPElectrostrictiveCouplingPolarDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (-2*(_q11[_qp]*_disp_x_grad[_qp](0) + _q12[_qp]*(_disp_y_grad[_qp](1) + _disp_z_grad[_qp](2))));
  }
  else if (_component == 1)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (-2*(_q11[_qp]*_disp_y_grad[_qp](1) + _q12[_qp]*(_disp_x_grad[_qp](0) + _disp_z_grad[_qp](2))));
  }
  else if (_component == 2)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (-2*(_q12[_qp]*(_disp_x_grad[_qp](0) + _disp_y_grad[_qp](1)) + _q11[_qp]*_disp_z_grad[_qp](2)));
  }
  else
    return 0.0;
}

Real
InducedPElectrostrictiveCouplingPolarDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _induced_polar_y_var)
    {
      return -_test[_i][_qp] * _phi[_j][_qp] * (-(_q44[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)))/2.);
    }
    else if (jvar == _induced_polar_z_var)
    {
      return -_test[_i][_qp] * _phi[_j][_qp] * (-(_q44[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))/2.);
    }
    /*else if (jvar == _disp_x_var)
    {
      return -_test[_i][_qp] * (0.0);
    }
    else if (jvar == _disp_y_var)
    {
      return -_test[_i][_qp] * (0.0);
    }
    else if (jvar == _disp_z_var)
    {
      return -_test[_i][_qp] * (0.0);
    }*/
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _induced_polar_x_var)
    {
      return -_test[_i][_qp] * (-(_q44[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)))/2.);
    }
    else if (jvar == _induced_polar_z_var)
    {
      return -_test[_i][_qp] * (-(_q44[_qp]*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))/2.);
    }
    /*else if (jvar == _disp_x_var)
    {
      return -_test[_i][_qp] * (0.0);
    }
    else if (jvar == _disp_y_var)
    {
      return -_test[_i][_qp] * (0.0);
    }
    else if (jvar == _disp_z_var)
    {
      return -_test[_i][_qp] * (0.0);
    }*/
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _induced_polar_x_var)
    {
      return -_test[_i][_qp] * (-(_q44[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))/2.);
    }
    else if (jvar == _induced_polar_y_var)
    {
      return -_test[_i][_qp] * (-(_q44[_qp]*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))/2.);
    }
    /*else if (jvar == _disp_x_var)
    {
      return -_test[_i][_qp] * (0.0);
    }
    else if (jvar == _disp_y_var)
    {
      return -_test[_i][_qp] * (0.0);
    }
    else if (jvar == _disp_z_var)
    {
      return -_test[_i][_qp] * (0.0);
    }*/
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
