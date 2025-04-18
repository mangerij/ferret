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

#include "UniaxialAFMSublattice.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", UniaxialAFMSublattice);

InputParameters UniaxialAFMSublattice::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for an uniaxial AFM sublattice");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredParam<unsigned int>("mag_sub", "An integer corresponding to the sublattice this Kernel acts on");
  params.addRequiredCoupledVar("mag1_x", "The x component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_y", "The y component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_z", "The z component of the constrained 1st sublattice magnetization vector");
  params.addCoupledVar("mag2_x", 0.0, "The x component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_y", 0.0, "The y component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_z", 0.0, "The z component of the constrained 2nd sublattice magnetization vector");
  return params;
}

UniaxialAFMSublattice::UniaxialAFMSublattice(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _mag_sub(getParam<unsigned int>("mag_sub")),
   _mag1_x_var(coupled("mag1_x")),
   _mag1_y_var(coupled("mag1_y")),
   _mag1_z_var(coupled("mag1_z")),
   _mag1_x(coupledValue("mag1_x")),
   _mag1_y(coupledValue("mag1_y")),
   _mag1_z(coupledValue("mag1_z")),
   _mag2_x_var(coupled("mag2_x")),
   _mag2_y_var(coupled("mag2_y")),
   _mag2_z_var(coupled("mag2_z")),
   _mag2_x(coupledValue("mag2_x")),
   _mag2_y(coupledValue("mag2_y")),
   _mag2_z(coupledValue("mag2_z")),
   _g0(getMaterialProperty<Real>("g0")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _H0(getMaterialProperty<Real>("H0")),
   _He(getMaterialProperty<Real>("He")),
   _Ha(getMaterialProperty<Real>("Ha"))
{
}


Real
UniaxialAFMSublattice::computeQpResidual()
{
  if (_mag_sub == 0)
  {
    if (_component == 0)
    {
      return _test[_i][_qp]*(_g0[_qp]*_Ms[_qp]*(_H0[_qp]*_mag1_y[_qp] + _He[_qp]*_mag1_z[_qp]*_mag2_y[_qp] - _He[_qp]*_mag1_y[_qp]*_mag2_z[_qp] + _Ha[_qp]*_mag1_y[_qp]*_mag1_z[_qp]*Utility::pow<2>(_Ms[_qp])));
    }
    else if (_component == 1)
    {
      return _test[_i][_qp]*(-(_g0[_qp]*_Ms[_qp]*(_H0[_qp]*_mag1_x[_qp] + _He[_qp]*_mag1_z[_qp]*_mag2_x[_qp] - _He[_qp]*_mag1_x[_qp]*_mag2_z[_qp] + _Ha[_qp]*_mag1_x[_qp]*_mag1_z[_qp]*Utility::pow<2>(_Ms[_qp]))));
    }
    else if (_component == 2)
    {
      return _test[_i][_qp]*(_g0[_qp]*_He[_qp]*(_mag1_y[_qp]*_mag2_x[_qp] - _mag1_x[_qp]*_mag2_y[_qp])*_Ms[_qp]);
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    if (_component == 0)
    {
      return _test[_i][_qp]*(_g0[_qp]*_Ms[_qp]*(_H0[_qp]*_mag2_y[_qp] - _He[_qp]*_mag1_z[_qp]*_mag2_y[_qp] + _He[_qp]*_mag1_y[_qp]*_mag2_z[_qp] + _Ha[_qp]*_mag2_y[_qp]*_mag2_z[_qp]*Utility::pow<2>(_Ms[_qp])));
    }
    else if (_component == 1)
    {
      return _test[_i][_qp]*(-(_g0[_qp]*_Ms[_qp]*(_H0[_qp]*_mag2_x[_qp] - _He[_qp]*_mag1_z[_qp]*_mag2_x[_qp] + _He[_qp]*_mag1_x[_qp]*_mag2_z[_qp] + _Ha[_qp]*_mag2_x[_qp]*_mag2_z[_qp]*Utility::pow<2>(_Ms[_qp]))));
    }
    else if (_component == 2)
    {
      return _test[_i][_qp]*(_g0[_qp]*_He[_qp]*(-(_mag1_y[_qp]*_mag2_x[_qp]) + _mag1_x[_qp]*_mag2_y[_qp])*_Ms[_qp]);
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
UniaxialAFMSublattice::computeQpJacobian()
{
//
//Note all on-diag jacobian entries are zero.
//
  if (_mag_sub == 0)
  {
    if (_component == 0)
    {
      return 0.0;
    }
    else if (_component == 1)
    {
      return 0.0;
    }
    else if (_component == 2)
    {
      return 0.0;
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    if (_component == 0)
    {
      return 0.0;
    }
    else if (_component == 1)
    {
      return 0.0;
    }
    else if (_component == 2)
    {
      return 0.0;
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
UniaxialAFMSublattice::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_mag_sub == 0)
  {
    if (_component == 0)
    {
      if (jvar == _mag1_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_Ms[_qp]*(_H0[_qp] - _He[_qp]*_mag2_z[_qp] + _Ha[_qp]*_mag1_z[_qp]*Utility::pow<2>(_Ms[_qp])));
      }
      else if (jvar == _mag1_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_Ms[_qp]*(_He[_qp]*_mag2_y[_qp] + _Ha[_qp]*_mag1_y[_qp]*Utility::pow<2>(_Ms[_qp])));
      }
      else if (jvar == _mag2_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag1_z[_qp]*_Ms[_qp]);
      }
      else if (jvar == _mag2_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag1_y[_qp]*_Ms[_qp]));
      }
      else
        return 0.0;
    }
    else if (_component == 1)
    {
      if (jvar == _mag1_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_Ms[_qp]*(_H0[_qp] - _He[_qp]*_mag2_z[_qp] + _Ha[_qp]*_mag1_z[_qp]*Utility::pow<2>(_Ms[_qp]))));
      }
      else if (jvar == _mag1_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_Ms[_qp]*(_He[_qp]*_mag2_x[_qp] + _Ha[_qp]*_mag1_x[_qp]*Utility::pow<2>(_Ms[_qp]))));
      }
      else if (jvar == _mag2_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag1_z[_qp]*_Ms[_qp]));
      }
      else if (jvar == _mag2_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag1_x[_qp]*_Ms[_qp]);
      }
      else
        return 0.0;
    }
    else if (_component == 2)
    {
      if (jvar == _mag1_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag2_y[_qp]*_Ms[_qp]));
      }
      else if (jvar == _mag1_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag2_x[_qp]*_Ms[_qp]);
      }
      else if (jvar == _mag2_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag1_y[_qp]*_Ms[_qp]);
      }
      else if (jvar == _mag2_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag1_x[_qp]*_Ms[_qp]));
      }
      else if (jvar == _mag2_z_var)
      {
        return 0.0;
      }
      else
        return 0.0;
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    if (_component == 0)
    {
      if (jvar == _mag2_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_Ms[_qp]*(_H0[_qp] - _He[_qp]*_mag1_z[_qp] + _Ha[_qp]*_mag2_z[_qp]*Utility::pow<2>(_Ms[_qp])));
      }
      else if (jvar == _mag2_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_Ms[_qp]*(_He[_qp]*_mag1_y[_qp] + _Ha[_qp]*_mag2_y[_qp]*Utility::pow<2>(_Ms[_qp])));
      }
      else if (jvar == _mag1_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag2_z[_qp]*_Ms[_qp]);
      }
      else if (jvar == _mag1_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag2_y[_qp]*_Ms[_qp]));
      }
      else
        return 0.0;
    }
    else if (_component == 1)
    {
      if (jvar == _mag2_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_Ms[_qp]*(_H0[_qp] - _He[_qp]*_mag1_z[_qp] + _Ha[_qp]*_mag2_z[_qp]*Utility::pow<2>(_Ms[_qp]))));
      }
      else if (jvar == _mag2_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_Ms[_qp]*(_He[_qp]*_mag1_x[_qp] + _Ha[_qp]*_mag2_x[_qp]*Utility::pow<2>(_Ms[_qp]))));
      }
      else if (jvar == _mag1_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag2_z[_qp]*_Ms[_qp]));
      }
      else if (jvar == _mag1_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_z_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag2_x[_qp]*_Ms[_qp]);
      }
      else
        return 0.0;
    }
    else if (_component == 2)
    {
      if (jvar == _mag2_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag1_y[_qp]*_Ms[_qp]));
      }
      else if (jvar == _mag2_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag1_x[_qp]*_Ms[_qp]);
      }
      else if (jvar == _mag1_x_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(_g0[_qp]*_He[_qp]*_mag2_y[_qp]*_Ms[_qp]);
      }
      else if (jvar == _mag1_y_var)
      {
        return _test[_i][_qp]*_phi[_j][_qp]*(-(_g0[_qp]*_He[_qp]*_mag2_x[_qp]*_Ms[_qp]));
      }
      else if (jvar == _mag1_z_var)
      {
        return 0.0;
      }
      else
        return 0.0;
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
