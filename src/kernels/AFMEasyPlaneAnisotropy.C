/*
   This file is part of FERRET, an add-on module for MOOSE
   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) af(1) later version.
   This program is distributed in the hope that it will be useful,
   but WITHOUT Af(1) WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret
**/

#include "AFMEasyPlaneAnisotropy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AFMEasyPlaneAnisotropy);

InputParameters AFMEasyPlaneAnisotropy::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredParam<unsigned int>("mag_sub", "An integer corresponding to the sublattice this Kernel acts on");
  params.addRequiredCoupledVar("mag1_x", "The x component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_y", "The y component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_z", "The z component of the constrained 1st sublattice magnetization vector");
  params.addCoupledVar("mag2_x", 0.0, "The x component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_y", 0.0, "The y component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_z", 0.0, "The z component of the constrained 2nd sublattice magnetization vector");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

AFMEasyPlaneAnisotropy::AFMEasyPlaneAnisotropy(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_sub(getParam<unsigned int>("mag_sub")),
  _mag1_x_var(coupled("mag1_x")),
  _mag1_y_var(coupled("mag1_y")),
  _mag1_z_var(coupled("mag1_z")),
  _mag2_x_var(coupled("mag2_x")),
  _mag2_y_var(coupled("mag2_y")),
  _mag2_z_var(coupled("mag2_z")),
  _mag1_x(coupledValue("mag1_x")),
  _mag1_y(coupledValue("mag1_y")),
  _mag1_z(coupledValue("mag1_z")),
  _mag2_x(coupledValue("mag2_x")),
  _mag2_y(coupledValue("mag2_y")),
  _mag2_z(coupledValue("mag2_z")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _K1(getMaterialProperty<Real>("K1")),
  _g0mu0Ms(getMaterialProperty<Real>("g0mu0Ms"))
{
}

Real
AFMEasyPlaneAnisotropy::computeQpResidual()
{
  if (_mag_sub == 0)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    if (_component == 0)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2))*(-(_mag1_z[_qp]*f(1)) + _mag1_y[_qp]*f(2) + _alpha[_qp]*(-(Utility::pow<2>(_mag1_y[_qp])*f(0)) + _mag1_x[_qp]*_mag1_y[_qp]*f(1) + _mag1_z[_qp]*(-(_mag1_z[_qp]*f(0)) + _mag1_x[_qp]*f(2))))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2))*(_alpha[_qp]*Utility::pow<2>(_mag1_z[_qp])*f(1) + _mag1_x[_qp]*(-(_alpha[_qp]*_mag1_y[_qp]*f(0)) + _alpha[_qp]*_mag1_x[_qp]*f(1) + f(2)) - _mag1_z[_qp]*(f(0) + _alpha[_qp]*_mag1_y[_qp]*f(2)))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2))*(_mag1_y[_qp]*(f(0) - _alpha[_qp]*_mag1_z[_qp]*f(1)) + _alpha[_qp]*Utility::pow<2>(_mag1_y[_qp])*f(2) + _mag1_x[_qp]*(-(_alpha[_qp]*_mag1_z[_qp]*f(0)) - f(1) + _alpha[_qp]*_mag1_x[_qp]*f(2)))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    if (_component == 0)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag2_x[_qp]*f(0) + _mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2))*(-(_mag2_z[_qp]*f(1)) + _mag2_y[_qp]*f(2) + _alpha[_qp]*(-(Utility::pow<2>(_mag2_y[_qp])*f(0)) + _mag2_x[_qp]*_mag2_y[_qp]*f(1) + _mag2_z[_qp]*(-(_mag2_z[_qp]*f(0)) + _mag2_x[_qp]*f(2))))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag2_x[_qp]*f(0) + _mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2))*(_alpha[_qp]*Utility::pow<2>(_mag2_z[_qp])*f(1) + _mag2_x[_qp]*(-(_alpha[_qp]*_mag2_y[_qp]*f(0)) + _alpha[_qp]*_mag2_x[_qp]*f(1) + f(2)) - _mag2_z[_qp]*(f(0) + _alpha[_qp]*_mag2_y[_qp]*f(2)))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag2_x[_qp]*f(0) + _mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2))*(_mag2_y[_qp]*(f(0) - _alpha[_qp]*_mag2_z[_qp]*f(1)) + _alpha[_qp]*Utility::pow<2>(_mag2_y[_qp])*f(2) + _mag2_x[_qp]*(-(_alpha[_qp]*_mag2_z[_qp]*f(0)) - f(1) + _alpha[_qp]*_mag2_x[_qp]*f(2)))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
AFMEasyPlaneAnisotropy::computeQpJacobian()
{
  if (_mag_sub == 0)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    if (_component == 0)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-(_mag1_z[_qp]*f(0)*f(1)) + _mag1_y[_qp]*f(0)*f(2) + _alpha[_qp]*(Utility::pow<2>(_mag1_y[_qp])*(-Utility::pow<2>(f(0)) + Utility::pow<2>(f(1))) + 2.0*_mag1_y[_qp]*f(1)*(_mag1_x[_qp]*f(0) + _mag1_z[_qp]*f(2)) + _mag1_z[_qp]*(-(_mag1_z[_qp]*Utility::pow<2>(f(0))) + 2.0*_mag1_x[_qp]*f(0)*f(2) + _mag1_z[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(f(1)*(_mag1_z[_qp]*f(0) - _mag1_x[_qp]*f(2)) + _alpha[_qp]*(Utility::pow<2>(_mag1_x[_qp])*(Utility::pow<2>(f(0)) - Utility::pow<2>(f(1))) + 2.0*_mag1_x[_qp]*f(0)*(_mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2)) + _mag1_z[_qp]*(-(_mag1_z[_qp]*Utility::pow<2>(f(1))) + 2.0*_mag1_y[_qp]*f(1)*f(2) + _mag1_z[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-(_mag1_y[_qp]*f(0)*f(2)) + _mag1_x[_qp]*f(1)*f(2) + _alpha[_qp]*(2.0*_mag1_x[_qp]*f(0)*(_mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2)) + Utility::pow<2>(_mag1_x[_qp])*(Utility::pow<2>(f(0)) - Utility::pow<2>(f(2))) + _mag1_y[_qp]*(_mag1_y[_qp]*Utility::pow<2>(f(1)) + 2.0*_mag1_z[_qp]*f(1)*f(2) - _mag1_y[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    if (_component == 0)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-(_mag2_z[_qp]*f(0)*f(1)) + _mag2_y[_qp]*f(0)*f(2) + _alpha[_qp]*(Utility::pow<2>(_mag2_y[_qp])*(-Utility::pow<2>(f(0)) + Utility::pow<2>(f(1))) + 2.0*_mag2_y[_qp]*f(1)*(_mag2_x[_qp]*f(0) + _mag2_z[_qp]*f(2)) + _mag2_z[_qp]*(-(_mag2_z[_qp]*Utility::pow<2>(f(0))) + 2.0*_mag2_x[_qp]*f(0)*f(2) + _mag2_z[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(f(1)*(_mag2_z[_qp]*f(0) - _mag2_x[_qp]*f(2)) + _alpha[_qp]*(Utility::pow<2>(_mag2_x[_qp])*(Utility::pow<2>(f(0)) - Utility::pow<2>(f(1))) + 2.0*_mag2_x[_qp]*f(0)*(_mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2)) + _mag2_z[_qp]*(-(_mag2_z[_qp]*Utility::pow<2>(f(1))) + 2.0*_mag2_y[_qp]*f(1)*f(2) + _mag2_z[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return -(-2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-(_mag2_y[_qp]*f(0)*f(2)) + _mag2_x[_qp]*f(1)*f(2) + _alpha[_qp]*(2.0*_mag2_x[_qp]*f(0)*(_mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2)) + Utility::pow<2>(_mag2_x[_qp])*(Utility::pow<2>(f(0)) - Utility::pow<2>(f(2))) + _mag2_y[_qp]*(_mag2_y[_qp]*Utility::pow<2>(f(1)) + 2.0*_mag2_z[_qp]*f(1)*f(2) - _mag2_y[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
AFMEasyPlaneAnisotropy::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_mag_sub == 0)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    if (_component == 0)
    {
      if (jvar == _mag1_y_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-((-2.0*_alpha[_qp]*_mag1_y[_qp]*f(0) + _alpha[_qp]*_mag1_x[_qp]*f(1) + f(2))*(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2))) + f(1)*(_mag1_z[_qp]*f(1) - _mag1_y[_qp]*f(2) + _alpha[_qp]*(Utility::pow<2>(_mag1_y[_qp])*f(0) - _mag1_x[_qp]*_mag1_y[_qp]*f(1) + _mag1_z[_qp]*(_mag1_z[_qp]*f(0) - _mag1_x[_qp]*f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_z_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag1_x[_qp]*f(0)*f(1) + _mag1_y[_qp]*Utility::pow<2>(f(1)) + 2.0*_mag1_z[_qp]*f(1)*f(2) - _mag1_y[_qp]*Utility::pow<2>(f(2)) + _alpha[_qp]*(-(Utility::pow<2>(_mag1_x[_qp])*f(0)*f(2)) + f(0)*(2.0*_mag1_y[_qp]*_mag1_z[_qp]*f(1) + Utility::pow<2>(_mag1_y[_qp])*f(2) + 3*Utility::pow<2>(_mag1_z[_qp])*f(2)) - 2.0*_mag1_x[_qp]*(_mag1_y[_qp]*f(1)*f(2) + _mag1_z[_qp]*(-Utility::pow<2>(f(0)) + Utility::pow<2>(f(2))))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_z_var)
      {
        return 0.0;
      }
      else
        return 0.0;
    }
    else if (_component == 1)
    {
      if (jvar == _mag1_x_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*((-(_alpha[_qp]*_mag1_y[_qp]*f(0)) + 2.0*_alpha[_qp]*_mag1_x[_qp]*f(1) + f(2))*(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2)) + f(0)*(_alpha[_qp]*Utility::pow<2>(_mag1_z[_qp])*f(1) + _mag1_x[_qp]*(-(_alpha[_qp]*_mag1_y[_qp]*f(0)) + _alpha[_qp]*_mag1_x[_qp]*f(1) + f(2)) - _mag1_z[_qp]*(f(0) + _alpha[_qp]*_mag1_y[_qp]*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_z_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-((f(0) - 2.0*_alpha[_qp]*_mag1_z[_qp]*f(1) + _alpha[_qp]*_mag1_y[_qp]*f(2))*(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2))) + f(2)*(_alpha[_qp]*Utility::pow<2>(_mag1_z[_qp])*f(1) + _mag1_x[_qp]*(-(_alpha[_qp]*_mag1_y[_qp]*f(0)) + _alpha[_qp]*_mag1_x[_qp]*f(1) + f(2)) - _mag1_z[_qp]*(f(0) + _alpha[_qp]*_mag1_y[_qp]*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_z_var)
      {
        return 0.0;
      }
      else
        return 0.0;
    }
    else if (_component == 2)
    {
      if (jvar == _mag1_x_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-((_alpha[_qp]*_mag1_z[_qp]*f(0) + f(1) - 2.0*_alpha[_qp]*_mag1_x[_qp]*f(2))*(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2))) + f(0)*(_mag1_y[_qp]*(f(0) - _alpha[_qp]*_mag1_z[_qp]*f(1)) + _alpha[_qp]*Utility::pow<2>(_mag1_y[_qp])*f(2) + _mag1_x[_qp]*(-(_alpha[_qp]*_mag1_z[_qp]*f(0)) - f(1) + _alpha[_qp]*_mag1_x[_qp]*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_y_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_alpha[_qp]*Utility::pow<2>(_mag1_x[_qp])*f(1)*f(2) + 3.0*_alpha[_qp]*Utility::pow<2>(_mag1_y[_qp])*f(1)*f(2) + _mag1_z[_qp]*(f(0) - _alpha[_qp]*_mag1_z[_qp]*f(1))*f(2) + _mag1_x[_qp]*(Utility::pow<2>(f(0)) - 2.0*_alpha[_qp]*_mag1_z[_qp]*f(0)*f(1) - Utility::pow<2>(f(1)) + 2.0*_alpha[_qp]*_mag1_y[_qp]*f(0)*f(2)) + 2.0*_mag1_y[_qp]*(f(0)*f(1) + _alpha[_qp]*_mag1_z[_qp]*(-Utility::pow<2>(f(1)) + Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_y_var)
      {
        return 0.0;
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
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    if (_component == 0)
    {
      if (jvar == _mag2_y_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-((-2.0*_alpha[_qp]*_mag2_y[_qp]*f(0) + _alpha[_qp]*_mag2_x[_qp]*f(1) + f(2))*(_mag2_x[_qp]*f(0) + _mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2))) + f(1)*(_mag2_z[_qp]*f(1) - _mag2_y[_qp]*f(2) + _alpha[_qp]*(Utility::pow<2>(_mag2_y[_qp])*f(0) - _mag2_x[_qp]*_mag2_y[_qp]*f(1) + _mag2_z[_qp]*(_mag2_z[_qp]*f(0) - _mag2_x[_qp]*f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_z_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_mag2_x[_qp]*f(0)*f(1) + _mag2_y[_qp]*Utility::pow<2>(f(1)) + 2.0*_mag2_z[_qp]*f(1)*f(2) - _mag2_y[_qp]*Utility::pow<2>(f(2)) + _alpha[_qp]*(-(Utility::pow<2>(_mag2_x[_qp])*f(0)*f(2)) + f(0)*(2.0*_mag2_y[_qp]*_mag2_z[_qp]*f(1) + Utility::pow<2>(_mag2_y[_qp])*f(2) + 3*Utility::pow<2>(_mag2_z[_qp])*f(2)) - 2.0*_mag2_x[_qp]*(_mag2_y[_qp]*f(1)*f(2) + _mag2_z[_qp]*(-Utility::pow<2>(f(0)) + Utility::pow<2>(f(2))))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_z_var)
      {
        return 0.0;
      }
      else
        return 0.0;
    }
    else if (_component == 1)
    {
      if (jvar == _mag2_x_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*((-(_alpha[_qp]*_mag2_y[_qp]*f(0)) + 2.0*_alpha[_qp]*_mag2_x[_qp]*f(1) + f(2))*(_mag2_x[_qp]*f(0) + _mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2)) + f(0)*(_alpha[_qp]*Utility::pow<2>(_mag2_z[_qp])*f(1) + _mag2_x[_qp]*(-(_alpha[_qp]*_mag2_y[_qp]*f(0)) + _alpha[_qp]*_mag2_x[_qp]*f(1) + f(2)) - _mag2_z[_qp]*(f(0) + _alpha[_qp]*_mag2_y[_qp]*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_z_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-((f(0) - 2.0*_alpha[_qp]*_mag2_z[_qp]*f(1) + _alpha[_qp]*_mag2_y[_qp]*f(2))*(_mag2_x[_qp]*f(0) + _mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2))) + f(2)*(_alpha[_qp]*Utility::pow<2>(_mag2_z[_qp])*f(1) + _mag2_x[_qp]*(-(_alpha[_qp]*_mag2_y[_qp]*f(0)) + _alpha[_qp]*_mag2_x[_qp]*f(1) + f(2)) - _mag2_z[_qp]*(f(0) + _alpha[_qp]*_mag2_y[_qp]*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_z_var)
      {
        return 0.0;
      }
      else
        return 0.0;
    }
    else if (_component == 2)
    {
      if (jvar == _mag2_x_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(-((_alpha[_qp]*_mag2_z[_qp]*f(0) + f(1) - 2.0*_alpha[_qp]*_mag2_x[_qp]*f(2))*(_mag2_x[_qp]*f(0) + _mag2_y[_qp]*f(1) + _mag2_z[_qp]*f(2))) + f(0)*(_mag2_y[_qp]*(f(0) - _alpha[_qp]*_mag2_z[_qp]*f(1)) + _alpha[_qp]*Utility::pow<2>(_mag2_y[_qp])*f(2) + _mag2_x[_qp]*(-(_alpha[_qp]*_mag2_z[_qp]*f(0)) - f(1) + _alpha[_qp]*_mag2_x[_qp]*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_y_var)
      {
        return -(2.0*_g0mu0Ms[_qp]*_K1[_qp]*(_alpha[_qp]*Utility::pow<2>(_mag2_x[_qp])*f(1)*f(2) + 3.0*_alpha[_qp]*Utility::pow<2>(_mag2_y[_qp])*f(1)*f(2) + _mag2_z[_qp]*(f(0) - _alpha[_qp]*_mag2_z[_qp]*f(1))*f(2) + _mag2_x[_qp]*(Utility::pow<2>(f(0)) - 2.0*_alpha[_qp]*_mag2_z[_qp]*f(0)*f(1) - Utility::pow<2>(f(1)) + 2.0*_alpha[_qp]*_mag2_y[_qp]*f(0)*f(2)) + 2.0*_mag2_y[_qp]*(f(0)*f(1) + _alpha[_qp]*_mag2_z[_qp]*(-Utility::pow<2>(f(1)) + Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_y_var)
      {
        return 0.0;
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
