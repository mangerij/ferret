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

#include "CorrectionFactorPolarBFO.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CorrectionFactorPolarBFO);

InputParameters CorrectionFactorPolarBFO::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a correction factor for the local free energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the tilt");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the tilt");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the tilt");
  return params;
}

CorrectionFactorPolarBFO::CorrectionFactorPolarBFO(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _f0(getMaterialProperty<Real>("f0")),
   _f1(getMaterialProperty<Real>("f1")),
   _c0(getMaterialProperty<Real>("c0"))
{
}

Real
CorrectionFactorPolarBFO::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_x[_qp]*
   (-3*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*(1 - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_x[_qp])*(-1 + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_y[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) - 3*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 + 2*Utility::pow<2>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_y[_qp])*(-1 + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_z[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<4>(_polar_y[_qp]) - 2*Utility::pow<2>(_polar_z[_qp]) + 2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) - 3*Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 - 3*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))));
  }
  else
    return 0.0;
}

Real
CorrectionFactorPolarBFO::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*
   (-15*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*(1 - 3*Utility::pow<2>(_polar_z[_qp])) + 6*Utility::pow<2>(_polar_x[_qp])*(-1 + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*
   (Utility::pow<4>(_polar_x[_qp]) - 15*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 + 6*Utility::pow<2>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])) + 6*Utility::pow<2>(_polar_y[_qp])*(-1 + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*
   (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])*(2 + 5*Utility::pow<2>(_polar_z[_qp])) + Utility::pow<2>(_polar_y[_qp])*(1 + 6*Utility::pow<2>(_polar_z[_qp])) + Utility::pow<2>(_polar_x[_qp])*(1 - 3*Utility::pow<2>(_polar_y[_qp]) + 6*Utility::pow<2>(_polar_z[_qp]))));
  }
  else
    return 0.0;
}

Real
CorrectionFactorPolarBFO::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_x[_qp]*_polar_y[_qp]*
   (1 + 2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_x[_qp]*_polar_z[_qp]*
   (1 + 2*Utility::pow<2>(_polar_x[_qp]) - 3*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp])));
    }
    else if (jvar == _antiphase_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_x[_qp]*(-3*Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp])) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_x[_qp]*
   (-3*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*(1 - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_x[_qp])*(-1 + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_y[_qp]*(Utility::pow<4>(_antiphase_A_x[_qp]) - 3*Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(1 + 2*Utility::pow<2>(_antiphase_A_y[_qp]) - 3*Utility::pow<2>(_antiphase_A_z[_qp])) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_x[_qp]*
   (-3*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*(1 - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_x[_qp])*(-1 + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_z[_qp]*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]) - 2*Utility::pow<2>(_antiphase_A_z[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) - 3*Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_x[_qp]*
   (-3*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*(1 - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_x[_qp])*(-1 + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
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
      return _test[_i][_qp] * _phi[_j][_qp] * (4*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_x[_qp]*_polar_y[_qp]*
   (1 + 2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_y[_qp]*_polar_z[_qp]*
   (1 - 3*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp])));
    }
    else if (jvar == _antiphase_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_x[_qp]*(-3*Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp])) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_y[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) - 3*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 + 2*Utility::pow<2>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_y[_qp])*(-1 + Utility::pow<2>(_polar_z[_qp]))));
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_y[_qp]*(Utility::pow<4>(_antiphase_A_x[_qp]) - 3*Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(1 + 2*Utility::pow<2>(_antiphase_A_y[_qp]) - 3*Utility::pow<2>(_antiphase_A_z[_qp])) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_y[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) - 3*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 + 2*Utility::pow<2>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_y[_qp])*(-1 + Utility::pow<2>(_polar_z[_qp]))));
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_z[_qp]*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]) - 2*Utility::pow<2>(_antiphase_A_z[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) - 3*Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_y[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) - 3*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 + 2*Utility::pow<2>(_polar_y[_qp]) - 3*Utility::pow<2>(_polar_z[_qp])) + 2*Utility::pow<2>(_polar_y[_qp])*(-1 + Utility::pow<2>(_polar_z[_qp]))));
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
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_x[_qp]*_polar_z[_qp]*
   (1 + 2*Utility::pow<2>(_polar_x[_qp]) - 3*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp])));
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_c0[_qp] + _f0[_qp])*_polar_y[_qp]*_polar_z[_qp]*
   (1 - 3*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp])));
    }
    else if (jvar == _antiphase_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_x[_qp]*(-3*Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp])) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_z[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<4>(_polar_y[_qp]) - 2*Utility::pow<2>(_polar_z[_qp]) + 2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) - 3*Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 - 3*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))));
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_y[_qp]*(Utility::pow<4>(_antiphase_A_x[_qp]) - 3*Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(1 + 2*Utility::pow<2>(_antiphase_A_y[_qp]) - 3*Utility::pow<2>(_antiphase_A_z[_qp])) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_z[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<4>(_polar_y[_qp]) - 2*Utility::pow<2>(_polar_z[_qp]) + 2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) - 3*Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 - 3*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))));
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*_antiphase_A_z[_qp]*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]) - 2*Utility::pow<2>(_antiphase_A_z[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) - 3*Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_z[_qp])))*_c0[_qp]*_polar_z[_qp]*
   (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<4>(_polar_y[_qp]) - 2*Utility::pow<2>(_polar_z[_qp]) + 2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) - 3*Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(1 - 3*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
