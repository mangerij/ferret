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

#include "AFMSublatticeAnisotropy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AFMSublatticeAnisotropy);

InputParameters AFMSublatticeAnisotropy::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

AFMSublatticeAnisotropy::AFMSublatticeAnisotropy(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _K1(getMaterialProperty<Real>("K1")),
  _g0(getMaterialProperty<Real>("g0")),
  _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
AFMSublatticeAnisotropy::computeQpResidual()
{
  if (_component == 0)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    return -(2.0*_g0[_qp]*_K1[_qp]*(_mag_x[_qp]*f(0) + _mag_y[_qp]*f(1) + _mag_z[_qp]*f(2))*(-(_mag_z[_qp]*f(1)) + _mag_y[_qp]*f(2) + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*f(0)) + _mag_x[_qp]*_mag_y[_qp]*f(1) + _mag_z[_qp]*(-(_mag_z[_qp]*f(0)) + _mag_x[_qp]*f(2))))*_test[_i][_qp])/
   ((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
  }
  else if (_component == 1)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    return -(2.0*_g0[_qp]*_K1[_qp]*(_mag_x[_qp]*f(0) + _mag_y[_qp]*f(1) + _mag_z[_qp]*f(2))*(-(_mag_z[_qp]*f(1)) + _mag_y[_qp]*f(2) + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*f(0)) + _mag_x[_qp]*_mag_y[_qp]*f(1) + _mag_z[_qp]*(-(_mag_z[_qp]*f(0)) + _mag_x[_qp]*f(2))))*_test[_i][_qp])/
   ((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
  }
  else if (_component == 2)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    return -(2.0*_g0[_qp]*_K1[_qp]*(_mag_x[_qp]*f(0) + _mag_y[_qp]*f(1) + _mag_z[_qp]*f(2))*(-(_mag_z[_qp]*f(1)) + _mag_y[_qp]*f(2) + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*f(0)) + _mag_x[_qp]*_mag_y[_qp]*f(1) + _mag_z[_qp]*(-(_mag_z[_qp]*f(0)) + _mag_x[_qp]*f(2))))*_test[_i][_qp])/
   ((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
  }
  else
    return 0.0;
}

Real
AFMSublatticeAnisotropy::computeQpJacobian()
{
  if (_component == 0)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    return -(2.0*_g0[_qp]*_K1[_qp]*(-(_mag_z[_qp]*f(0)*f(1)) + _mag_y[_qp]*f(0)*f(2) + _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_y[_qp])*(-Utility::pow<2>(f(0)) + Utility::pow<2>(f(1))) + 2.0*_mag_y[_qp]*f(1)*(_mag_x[_qp]*f(0) + _mag_z[_qp]*f(2)) + 
          _mag_z[_qp]*(-(_mag_z[_qp]*Utility::pow<2>(f(0))) + 2.0*_mag_x[_qp]*f(0)*f(2) + _mag_z[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
  }
  else if (_component == 1)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    return -(2.0*_g0[_qp]*_K1[_qp]*(-(_mag_z[_qp]*Utility::pow<2>(f(1))) + _mag_x[_qp]*f(0)*f(2) + 2*_mag_y[_qp]*f(1)*f(2) + _mag_z[_qp]*Utility::pow<2>(f(2)) + 
       _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_x[_qp])*f(0)*f(1) - f(0)*(3*Utility::pow<2>(_mag_y[_qp])*f(1) + Utility::pow<2>(_mag_z[_qp])*f(1) + 2*_mag_y[_qp]*_mag_z[_qp]*f(2)) + 
          2.0*_mag_x[_qp]*(-(_mag_y[_qp]*Utility::pow<2>(f(0))) + _mag_y[_qp]*Utility::pow<2>(f(1)) + _mag_z[_qp]*f(1)*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
  }
  else if (_component == 2)
  {
    RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
    RealVectorValue f = w/std::sqrt(w*w);
    return -(2.0*_g0[_qp]*_K1[_qp]*(-((2*_alpha[_qp]*_Ms[_qp]*_mag_z[_qp]*f(0) + f(1))*(_mag_x[_qp]*f(0) + _mag_y[_qp]*f(1))) - 2.0*_mag_z[_qp]*f(1)*f(2) + 
       _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_x[_qp])*f(0) - (Utility::pow<2>(_mag_y[_qp]) + 3*Utility::pow<2>(_mag_z[_qp]))*f(0) + 2.0*_mag_x[_qp]*_mag_y[_qp]*f(1))*f(2) + (_mag_y[_qp] + 2.0*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp])*Utility::pow<2>(f(2)))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
  }
  else
    return 0.0;
}

Real
AFMSublatticeAnisotropy::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
      RealVectorValue f = w/std::sqrt(w*w);
      return -(2.0*_g0[_qp]*_K1[_qp]*(-(_mag_z[_qp]*Utility::pow<2>(f(1))) + _mag_x[_qp]*f(0)*f(2) + 2*_mag_y[_qp]*f(1)*f(2) + _mag_z[_qp]*Utility::pow<2>(f(2)) + 
       _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_x[_qp])*f(0)*f(1) - f(0)*(3.0*Utility::pow<2>(_mag_y[_qp])*f(1) + Utility::pow<2>(_mag_z[_qp])*f(1) + 2.0*_mag_y[_qp]*_mag_z[_qp]*f(2)) + 
          2.0*_mag_x[_qp]*(-(_mag_y[_qp]*Utility::pow<2>(f(0))) + _mag_y[_qp]*Utility::pow<2>(f(1)) + _mag_z[_qp]*f(1)*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
      RealVectorValue f = w/std::sqrt(w*w);
      return -(2.0*_g0[_qp]*_K1[_qp]*(-((2*_alpha[_qp]*_Ms[_qp]*_mag_z[_qp]*f(0) + f(1))*(_mag_x[_qp]*f(0) + _mag_y[_qp]*f(1))) - 2.0*_mag_z[_qp]*f(1)*f(2) + 
       _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_x[_qp])*f(0) - (Utility::pow<2>(_mag_y[_qp]) + 3*Utility::pow<2>(_mag_z[_qp]))*f(0) + 2.0*_mag_x[_qp]*_mag_y[_qp]*f(1))*f(2) + (_mag_y[_qp] + 2.0*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp])*Utility::pow<2>(f(2)))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
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
      RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
      RealVectorValue f = w/std::sqrt(w*w);
      return -(2.0*_g0[_qp]*_K1[_qp]*(-(_mag_z[_qp]*f(0)*f(1)) + _mag_y[_qp]*f(0)*f(2) + _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_y[_qp])*(-Utility::pow<2>(f(0)) + Utility::pow<2>(f(1))) + 2.0*_mag_y[_qp]*f(1)*(_mag_x[_qp]*f(0) + _mag_z[_qp]*f(2)) + 
          _mag_z[_qp]*(-(_mag_z[_qp]*Utility::pow<2>(f(0))) + 2.0*_mag_x[_qp]*f(0)*f(2) + _mag_z[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
      RealVectorValue f = w/std::sqrt(w*w);
      return -(2.0*_g0[_qp]*_K1[_qp]*(-((2.0*_alpha[_qp]*_Ms[_qp]*_mag_z[_qp]*f(0) + f(1))*(_mag_x[_qp]*f(0) + _mag_y[_qp]*f(1))) - 2*_mag_z[_qp]*f(1)*f(2) + 
       _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_x[_qp])*f(0) - (Utility::pow<2>(_mag_y[_qp]) + 3.0*Utility::pow<2>(_mag_z[_qp]))*f(0) + 2*_mag_x[_qp]*_mag_y[_qp]*f(1))*f(2) + (_mag_y[_qp] + 2.0*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp])*Utility::pow<2>(f(2)))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
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
      RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
      RealVectorValue f = w/std::sqrt(w*w);
      return -(2.0*_g0[_qp]*_K1[_qp]*(-(_mag_z[_qp]*f(0)*f(1)) + _mag_y[_qp]*f(0)*f(2) + _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_y[_qp])*(-Utility::pow<2>(f(0)) + Utility::pow<2>(f(1))) + 2*_mag_y[_qp]*f(1)*(_mag_x[_qp]*f(0) + _mag_z[_qp]*f(2)) + 
          _mag_z[_qp]*(-(_mag_z[_qp]*Utility::pow<2>(f(0))) + 2*_mag_x[_qp]*f(0)*f(2) + _mag_z[_qp]*Utility::pow<2>(f(2)))))*_phi[_j][_qp]*_test[_i][_qp])/((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
      RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
      RealVectorValue f = w/std::sqrt(w*w);
      return -(2.0*_g0[_qp]*_K1[_qp]*(-(_mag_z[_qp]*Utility::pow<2>(f(1))) + _mag_x[_qp]*f(0)*f(2) + 2*_mag_y[_qp]*f(1)*f(2) + _mag_z[_qp]*Utility::pow<2>(f(2)) + 
       _alpha[_qp]*_Ms[_qp]*(Utility::pow<2>(_mag_x[_qp])*f(0)*f(1) - f(0)*(3*Utility::pow<2>(_mag_y[_qp])*f(1) + Utility::pow<2>(_mag_z[_qp])*f(1) + 2.0*_mag_y[_qp]*_mag_z[_qp]*f(2)) + 
          2.0*_mag_x[_qp]*(-(_mag_y[_qp]*Utility::pow<2>(f(0))) + _mag_y[_qp]*Utility::pow<2>(f(1)) + _mag_z[_qp]*f(1)*f(2))))*_phi[_j][_qp]*_test[_i][_qp])/((1.0 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<3>(_Ms[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
