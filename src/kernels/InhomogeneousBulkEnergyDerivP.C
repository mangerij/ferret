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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "InhomogeneousBulkEnergyDerivP.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", InhomogeneousBulkEnergyDerivP);

InputParameters InhomogeneousBulkEnergyDerivP::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("x", "The concentration");
  params.addParam<Real>("alpha01", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha011", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha0111", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha012", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha0112", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha0123", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1111", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1112", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1122", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1123", "The coefficients of the Landau expansion");
  params.addParam<Real>("b1", "The coupling coefficients to Sr concentration");
  params.addParam<Real>("b2", "The coupling coefficients to Sr concentration");
  params.addParam<Real>("b3", "The coupling coefficients to Sr concentration");
  params.addParam<Real>("T", "The temperature");
  return params;
}

InhomogeneousBulkEnergyDerivP::InhomogeneousBulkEnergyDerivP(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _x(coupledValue("x")),
   _alpha01(getParam<Real>("alpha01")),
   _alpha011(getParam<Real>("alpha011")),
   _alpha0111(getParam<Real>("alpha0111")),
   _alpha012(getParam<Real>("alpha012")),
   _alpha0112(getParam<Real>("alpha0112")),
   _alpha0123(getParam<Real>("alpha0123")),
   _alpha1111(getParam<Real>("alpha1111")),
   _alpha1112(getParam<Real>("alpha1112")),
   _alpha1122(getParam<Real>("alpha1122")),
   _alpha1123(getParam<Real>("alpha1123")),
   _b1(getParam<Real>("b1")),
   _b2(getParam<Real>("b2")),
   _b3(getParam<Real>("b3")),
   _T(getParam<Real>("T"))
{
}

Real
InhomogeneousBulkEnergyDerivP::computeQpResidual()
{
  if (_component == 0)
  {
      return _test[_i][_qp]* (2*_polar_x[_qp]*(4*_alpha1111*Utility::pow<6>(_polar_x[_qp]) + _alpha1123*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*(2*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])) + 2*_alpha1122*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp])) +
     _alpha1112*(Utility::pow<6>(_polar_y[_qp]) + Utility::pow<6>(_polar_z[_qp]) + 3*Utility::pow<4>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))) + _alpha01*(-1.5727418732060157 + 1/(-0.5 + 0.5*std::exp(2.0*(160./_T))) - 1.*_b1*_x[_qp]) -
     2*_alpha011*Utility::pow<2>(_polar_x[_qp])*(-1 + _b2*_x[_qp]) - _alpha012*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*(-1 + _b2*_x[_qp]) + _alpha0123*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*(1 - _b3*_x[_qp]) +
     _alpha0112*(Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp]) + 2*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*(1 - _b3*_x[_qp]) - 3*_alpha0111*Utility::pow<4>(_polar_x[_qp])*(-1 + _b3*_x[_qp])));
     }
  else if (_component == 1)
  {
    return _test[_i][_qp]*(2*_polar_y[_qp]*(4*_alpha1111*Utility::pow<6>(_polar_y[_qp]) + _alpha1123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])) + 2*_alpha1122*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_z[_qp])) +
     _alpha1112*(Utility::pow<6>(_polar_x[_qp]) + Utility::pow<6>(_polar_z[_qp]) + 3*Utility::pow<4>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))) + _alpha01*(-1.5727418732060157 + 1/(-0.5 + 0.5*std::exp(2.0*(160./_T))) - 1.*_b1*_x[_qp]) -
     2*_alpha011*Utility::pow<2>(_polar_y[_qp])*(-1 + _b2*_x[_qp]) - _alpha012*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*(-1 + _b2*_x[_qp]) + _alpha0123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*(1 - _b3*_x[_qp]) -
     3*_alpha0111*Utility::pow<4>(_polar_y[_qp])*(-1 + _b3*_x[_qp]) - _alpha0112*(Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]))*(-1 + _b3*_x[_qp])));
     }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_polar_z[_qp]*(2*_alpha1122*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]) + 4*_alpha1111*Utility::pow<6>(_polar_z[_qp]) + _alpha1123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp])) +
     _alpha1112*(Utility::pow<6>(_polar_x[_qp]) + Utility::pow<6>(_polar_y[_qp]) + 3*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<4>(_polar_z[_qp])) + _alpha01*(-1.5727418732060157 + 1/(-0.5 + 0.5*std::exp(2.0*(160./_T))) - 1.*_b1*_x[_qp]) -
     _alpha012*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*(-1 + _b2*_x[_qp]) - 2*_alpha011*Utility::pow<2>(_polar_z[_qp])*(-1 + _b2*_x[_qp]) + _alpha0123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*(1 - _b3*_x[_qp]) -
     3*_alpha0111*Utility::pow<4>(_polar_z[_qp])*(-1 + _b3*_x[_qp]) - _alpha0112*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + 2*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*(-1 + _b3*_x[_qp])));
     }
  else
    return 0.0;
}

Real
InhomogeneousBulkEnergyDerivP::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * 2*(28*_alpha1111*Utility::pow<6>(_polar_x[_qp]) + _alpha1123*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*(6*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])) +
     6*_alpha1122*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp])) + _alpha1112*(Utility::pow<6>(_polar_y[_qp]) + Utility::pow<6>(_polar_z[_qp]) + 15*Utility::pow<4>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))) +
     _alpha01*(-1.5727418732060157 + 1/(-0.5 + 0.5*std::exp(2.0*(160./_T))) - 1.*_b1*_x[_qp]) - 6*_alpha011*Utility::pow<2>(_polar_x[_qp])*(-1 + _b2*_x[_qp]) -
     _alpha012*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*(-1 + _b2*_x[_qp]) + _alpha0123*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*(1 - _b3*_x[_qp]) +
     _alpha0112*(Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp]) + 6*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*(1 - _b3*_x[_qp]) - 15*_alpha0111*Utility::pow<4>(_polar_x[_qp])*(-1 + _b3*_x[_qp]));
     }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] *2*(28*_alpha1111*Utility::pow<6>(_polar_y[_qp]) + _alpha1123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_polar_x[_qp]) + 6*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])) +
     6*_alpha1122*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_z[_qp])) + _alpha1112*(Utility::pow<6>(_polar_x[_qp]) + Utility::pow<6>(_polar_z[_qp]) + 15*Utility::pow<4>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))) +
     _alpha01*(-1.5727418732060157 + 1/(-0.5 + 0.5*std::exp(2.0*(160./_T))) - 1.*_b1*_x[_qp]) - 6*_alpha011*Utility::pow<2>(_polar_y[_qp])*(-1 + _b2*_x[_qp]) -
     _alpha012*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*(-1 + _b2*_x[_qp]) + _alpha0123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*(1 - _b3*_x[_qp]) +
     _alpha0112*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_z[_qp]) + 6*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])))*(1 - _b3*_x[_qp]) - 15*_alpha0111*Utility::pow<4>(_polar_y[_qp])*(-1 + _b3*_x[_qp]));
    }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] *2*(6*_alpha1122*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]) + 28*_alpha1111*Utility::pow<6>(_polar_z[_qp]) +
     _alpha1123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + 6*Utility::pow<2>(_polar_z[_qp])) +
     _alpha1112*(Utility::pow<6>(_polar_x[_qp]) + Utility::pow<6>(_polar_y[_qp]) + 15*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<4>(_polar_z[_qp])) + _alpha01*(-1.5727418732060157 + 1/(-0.5 + 0.5*std::exp(2.0*(160./_T))) - 1.*_b1*_x[_qp]) -
     _alpha012*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*(-1 + _b2*_x[_qp]) - 6*_alpha011*Utility::pow<2>(_polar_z[_qp])*(-1 + _b2*_x[_qp]) + _alpha0123*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*(1 - _b3*_x[_qp]) +
     _alpha0112*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + 6*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*(1 - _b3*_x[_qp]) - 15*_alpha0111*Utility::pow<4>(_polar_z[_qp])*(-1 + _b3*_x[_qp]));
     }
  else
    return 0.0;
}

Real
InhomogeneousBulkEnergyDerivP::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *4*_polar_x[_qp]*_polar_y[_qp]*(_alpha012 + 2*_alpha0112*Utility::pow<2>(_polar_x[_qp]) + 3*_alpha1112*Utility::pow<4>(_polar_x[_qp]) + 3*_alpha1112*Utility::pow<4>(_polar_y[_qp]) + (_alpha0123 + 2*_alpha1123*Utility::pow<2>(_polar_x[_qp]))*Utility::pow<2>(_polar_z[_qp]) +
     _alpha1123*Utility::pow<4>(_polar_z[_qp]) - (_alpha012*_b2 + 2*_alpha0112*_b3*Utility::pow<2>(_polar_x[_qp]) + _alpha0123*_b3*Utility::pow<2>(_polar_z[_qp]))*_x[_qp] +
     2*Utility::pow<2>(_polar_y[_qp])*(_alpha0112 + 2*_alpha1122*Utility::pow<2>(_polar_x[_qp]) + _alpha1123*Utility::pow<2>(_polar_z[_qp]) - _alpha0112*_b3*_x[_qp]));
     }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * 4*_polar_x[_qp]*_polar_z[_qp]*(_alpha012 + 2*_alpha0112*Utility::pow<2>(_polar_x[_qp]) + 3*_alpha1112*Utility::pow<4>(_polar_x[_qp]) + (_alpha0123 + 2*_alpha1123*Utility::pow<2>(_polar_x[_qp]))*Utility::pow<2>(_polar_y[_qp]) + _alpha1123*Utility::pow<4>(_polar_y[_qp]) +
     3*_alpha1112*Utility::pow<4>(_polar_z[_qp]) - (_alpha012*_b2 + 2*_alpha0112*_b3*Utility::pow<2>(_polar_x[_qp]) + _alpha0123*_b3*Utility::pow<2>(_polar_y[_qp]))*_x[_qp] +
     2*Utility::pow<2>(_polar_z[_qp])*(_alpha0112 + 2*_alpha1122*Utility::pow<2>(_polar_x[_qp]) + _alpha1123*Utility::pow<2>(_polar_y[_qp]) - _alpha0112*_b3*_x[_qp]));
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
      return _test[_i][_qp] * _phi[_j][_qp]*4*_polar_x[_qp]*_polar_y[_qp]*(_alpha012 + 3*_alpha1112*Utility::pow<4>(_polar_x[_qp]) + 2*_alpha0112*Utility::pow<2>(_polar_y[_qp]) + 3*_alpha1112*Utility::pow<4>(_polar_y[_qp]) + (_alpha0123 + 2*_alpha1123*Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]) +
     _alpha1123*Utility::pow<4>(_polar_z[_qp]) - (_alpha012*_b2 + 2*_alpha0112*_b3*Utility::pow<2>(_polar_y[_qp]) + _alpha0123*_b3*Utility::pow<2>(_polar_z[_qp]))*_x[_qp] +
     2*Utility::pow<2>(_polar_x[_qp])*(_alpha0112 + 2*_alpha1122*Utility::pow<2>(_polar_y[_qp]) + _alpha1123*Utility::pow<2>(_polar_z[_qp]) - _alpha0112*_b3*_x[_qp]));
     }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * 4*_polar_y[_qp]*_polar_z[_qp]*(_alpha012 + _alpha0123*Utility::pow<2>(_polar_x[_qp]) + _alpha1123*Utility::pow<4>(_polar_x[_qp]) + 2*(_alpha0112 + _alpha1123*Utility::pow<2>(_polar_x[_qp]))*Utility::pow<2>(_polar_y[_qp]) + 3*_alpha1112*Utility::pow<4>(_polar_y[_qp]) +
     3*_alpha1112*Utility::pow<4>(_polar_z[_qp]) - (_alpha012*_b2 + _alpha0123*_b3*Utility::pow<2>(_polar_x[_qp]) + 2*_alpha0112*_b3*Utility::pow<2>(_polar_y[_qp]))*_x[_qp] +
     2*Utility::pow<2>(_polar_z[_qp])*(_alpha0112 + _alpha1123*Utility::pow<2>(_polar_x[_qp]) + 2*_alpha1122*Utility::pow<2>(_polar_y[_qp]) - _alpha0112*_b3*_x[_qp]));
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
      return _test[_i][_qp] * _phi[_j][_qp] *4*_polar_x[_qp]*_polar_z[_qp]*(_alpha012 + 3*_alpha1112*Utility::pow<4>(_polar_x[_qp]) + _alpha0123*Utility::pow<2>(_polar_y[_qp]) + _alpha1123*Utility::pow<4>(_polar_y[_qp]) + 2*(_alpha0112 + _alpha1123*Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]) +
     3*_alpha1112*Utility::pow<4>(_polar_z[_qp]) - (_alpha012*_b2 + _alpha0123*_b3*Utility::pow<2>(_polar_y[_qp]) + 2*_alpha0112*_b3*Utility::pow<2>(_polar_z[_qp]))*_x[_qp] +
     2*Utility::pow<2>(_polar_x[_qp])*(_alpha0112 + _alpha1123*Utility::pow<2>(_polar_y[_qp]) + 2*_alpha1122*Utility::pow<2>(_polar_z[_qp]) - _alpha0112*_b3*_x[_qp]));
     }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *4*_polar_y[_qp]*_polar_z[_qp]*(_alpha012 + _alpha0123*Utility::pow<2>(_polar_x[_qp]) + _alpha1123*Utility::pow<4>(_polar_x[_qp]) + 3*_alpha1112*Utility::pow<4>(_polar_y[_qp]) + 2*(_alpha0112 + _alpha1123*Utility::pow<2>(_polar_x[_qp]))*Utility::pow<2>(_polar_z[_qp]) +
     3*_alpha1112*Utility::pow<4>(_polar_z[_qp]) - (_alpha012*_b2 + _alpha0123*_b3*Utility::pow<2>(_polar_x[_qp]) + 2*_alpha0112*_b3*Utility::pow<2>(_polar_z[_qp]))*_x[_qp] +
     2*Utility::pow<2>(_polar_y[_qp])*(_alpha0112 + _alpha1123*Utility::pow<2>(_polar_x[_qp]) + 2*_alpha1122*Utility::pow<2>(_polar_z[_qp]) - _alpha0112*_b3*_x[_qp]));
     }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
