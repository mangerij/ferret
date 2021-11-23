/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the ter___Ms of the GNU General Public License as published by
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

#include "SublatticeAFMCoupling.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", SublatticeAFMCoupling);

InputParameters SublatticeAFMCoupling::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("mag_x", "The x component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag2_x", "The x component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag2_y", "The y component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag2_z", "The z component of a constrained sublattice magnetic vector");
  return params;
}

SublatticeAFMCoupling::SublatticeAFMCoupling(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _mag2_x_var(coupled("mag2_x")),
  _mag2_y_var(coupled("mag2_y")),
  _mag2_z_var(coupled("mag2_z")),
  _mag2_x(coupledValue("mag2_x")),
  _mag2_y(coupledValue("mag2_y")),
  _mag2_z(coupledValue("mag2_z")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _g0(getMaterialProperty<Real>("g0")),
  _La(getMaterialProperty<Real>("La")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _mu0(getMaterialProperty<Real>("mu0"))
{
}

Real
SublatticeAFMCoupling::computeQpResidual()
{
  if (_component == 0)
  {
   return (_g0[_qp]*_La[_qp]*(_mag2_y[_qp]*_mag_z[_qp] + _mag2_z[_qp]*(-_mag_y[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp]) + _alpha[_qp]*_Ms[_qp]*(_mag2_y[_qp]*_mag_x[_qp]*_mag_y[_qp] - _mag2_x[_qp]*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
   return (_g0[_qp]*_La[_qp]*(-(_mag2_x[_qp]*_mag_z[_qp]) + _mag2_z[_qp]*(_mag_x[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]) - _alpha[_qp]*_Ms[_qp]*(-(_mag2_x[_qp]*_mag_x[_qp]*_mag_y[_qp]) + _mag2_y[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 2)
  {
   return (_g0[_qp]*_La[_qp]*(_mag2_x[_qp]*_mag_y[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag2_z[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))) + _mag2_x[_qp]*_mag_x[_qp]*_mag_z[_qp]) + _mag2_y[_qp]*(-_mag_x[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else
    return 0.0;
}

Real
SublatticeAFMCoupling::computeQpJacobian()
{
  if (_component == 0)
  {
  return (_alpha[_qp]*_g0[_qp]*_La[_qp]*(_mag2_y[_qp]*_mag_y[_qp] + _mag2_z[_qp]*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 1)
  {
  return (_alpha[_qp]*_g0[_qp]*_La[_qp]*(_mag2_x[_qp]*_mag_x[_qp] + _mag2_z[_qp]*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 2)
  {
  return (_alpha[_qp]*_g0[_qp]*_La[_qp]*(_mag2_x[_qp]*_mag_x[_qp] + _mag2_y[_qp]*_mag_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else
    return 0.0;
}

Real
SublatticeAFMCoupling::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
    return (_g0[_qp]*_La[_qp]*(-_mag2_z[_qp] + _alpha[_qp]*_Ms[_qp]*(_mag2_y[_qp]*_mag_x[_qp] - 2*_mag2_x[_qp]*_mag_y[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
    return (_g0[_qp]*_La[_qp]*(_mag2_y[_qp] + _alpha[_qp]*_Ms[_qp]*(_mag2_z[_qp]*_mag_x[_qp] - 2*_mag2_x[_qp]*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag2_x_var)
    {
    return -((_alpha[_qp]*_g0[_qp]*_La[_qp]*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1 + Utility::pow<2>(_alpha[_qp])));
    }
    else if (jvar == _mag2_y_var)
    {
    return (_g0[_qp]*_La[_qp]*(_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp] + _mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag2_z_var)
    {
    return (_g0[_qp]*_La[_qp]*(-_mag_y[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _mag_x_var)
    {
    return (_g0[_qp]*_La[_qp]*(_mag2_z[_qp] + _alpha[_qp]*_Ms[_qp]*(-2*_mag2_y[_qp]*_mag_x[_qp] + _mag2_x[_qp]*_mag_y[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]); 
    }
    else if (jvar == _mag_z_var)
    {
    return (_g0[_qp]*_La[_qp]*(-_mag2_x[_qp] + _alpha[_qp]*_Ms[_qp]*(_mag2_z[_qp]*_mag_y[_qp] - 2*_mag2_y[_qp]*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag2_x_var)
    {
    return (_g0[_qp]*_La[_qp]*(_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp] - _mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag2_y_var)
    {
    return -((_alpha[_qp]*_g0[_qp]*_La[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1 + Utility::pow<2>(_alpha[_qp])));
    }
    else if (jvar == _mag2_z_var)
    {
    return (_g0[_qp]*_La[_qp]*(_mag_x[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _mag_x_var)
    {
    return -((_g0[_qp]*_La[_qp]*(_mag2_y[_qp] + _alpha[_qp]*_Ms[_qp]*(2*_mag2_z[_qp]*_mag_x[_qp] - _mag2_x[_qp]*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
    return (_g0[_qp]*_La[_qp]*(_mag2_x[_qp] + _alpha[_qp]*_Ms[_qp]*(-2*_mag2_z[_qp]*_mag_y[_qp] + _mag2_y[_qp]*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag2_x_var)
    {
    return (_g0[_qp]*_La[_qp]*(_mag_y[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag2_y_var)
    {
    return (_g0[_qp]*_La[_qp]*(-_mag_x[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag2_z_var)
    {
    return -((_alpha[_qp]*_g0[_qp]*_La[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/(1 + Utility::pow<2>(_alpha[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
