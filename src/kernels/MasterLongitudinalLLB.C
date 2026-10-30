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

#include "MasterLongitudinalLLB.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MasterLongitudinalLLB);


InputParameters MasterLongitudinalLLB::validParams()

{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addParam<Real>("g0", 1.0, "electron gyromagnetic factor");
  params.addParam<Real>("fD", 1.0, "the magnitude of the normalized magnetization");
  params.addParam<MooseEnum>("form", MooseEnum("penalty llb", "penalty"), "penalty: the T = 0 norm penalty alpha_long |m|^2(|m|^2 - fD) m (default); llb: the Landau-Lifshitz-Bloch longitudinal field H_L = (1 - |m|^2/me^2) m/(2 chi_par) below Tc and -(1 + 3 Tc |m|^2/(5 (T - Tc))) m/chi_par above, acting through -alpha_par H_L (material properties me, chi_par [per field unit], alpha_par)");
  params.addCoupledVar("temperature", 0.0, "The temperature (form = llb)");
  params.addParam<Real>("Tc", 0.0, "The Curie temperature (form = llb)");
  params.addParam<Real>("Hscale", 1.0, "scaling factor for effective fields (form = llb)");
  return params;
}

MasterLongitudinalLLB::MasterLongitudinalLLB(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _g0(getParam<Real>("g0")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _alpha_long(getOptionalMaterialProperty<Real>("alpha_long")),
  _fD(getParam<Real>("fD")),
  _llb(getParam<MooseEnum>("form") == "llb"),
  _temperature(coupledValue("temperature")),
  _Tc(getParam<Real>("Tc")),
  _Hscale(getParam<Real>("Hscale")),
  _me(getOptionalMaterialProperty<Real>("me")),
  _chi_par(getOptionalMaterialProperty<Real>("chi_par")),
  _alpha_par(getOptionalMaterialProperty<Real>("alpha_par"))
{
}

Real
MasterLongitudinalLLB::computeQpResidual()
{
  if (_llb)
  {
    if (!_me || !_chi_par || !_alpha_par)
      mooseError("form = llb needs the material properties me, chi_par and alpha_par");
    if (_temperature[_qp] < _Tc)
    {
      if (_component == 0)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((1.0/2.0)*(_test[_i][_qp]*Utility::pow<2>(_me[_qp])*_mag_x[_qp] - _test[_i][_qp]*Utility::pow<3>(_mag_x[_qp]) - _test[_i][_qp]*_mag_x[_qp]*Utility::pow<2>(_mag_y[_qp]) - _test[_i][_qp]*_mag_x[_qp]*Utility::pow<2>(_mag_z[_qp]))/Utility::pow<2>(_me[_qp]));
      }
      else if (_component == 1)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((1.0/2.0)*(_test[_i][_qp]*Utility::pow<2>(_me[_qp])*_mag_y[_qp] - _test[_i][_qp]*Utility::pow<2>(_mag_x[_qp])*_mag_y[_qp] - _test[_i][_qp]*Utility::pow<3>(_mag_y[_qp]) - _test[_i][_qp]*_mag_y[_qp]*Utility::pow<2>(_mag_z[_qp]))/Utility::pow<2>(_me[_qp]));
      }
      else if (_component == 2)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((1.0/2.0)*(_test[_i][_qp]*Utility::pow<2>(_me[_qp])*_mag_z[_qp] - _test[_i][_qp]*Utility::pow<2>(_mag_x[_qp])*_mag_z[_qp] - _test[_i][_qp]*Utility::pow<2>(_mag_y[_qp])*_mag_z[_qp] - _test[_i][_qp]*Utility::pow<3>(_mag_z[_qp]))/Utility::pow<2>(_me[_qp]));
      }
      else
        return 0.0;
    }
    else
    {
      if (_component == 0)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((3*_test[_i][_qp]*_Tc*Utility::pow<3>(_mag_x[_qp]) + 3*_test[_i][_qp]*_Tc*_mag_x[_qp]*Utility::pow<2>(_mag_y[_qp]) + 3*_test[_i][_qp]*_Tc*_mag_x[_qp]*Utility::pow<2>(_mag_z[_qp]) - 5*_test[_i][_qp]*_Tc*_mag_x[_qp] + 5*_test[_i][_qp]*_temperature[_qp]*_mag_x[_qp])/(5*_Tc - 5*_temperature[_qp]));
      }
      else if (_component == 1)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((3*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_x[_qp])*_mag_y[_qp] + 3*_test[_i][_qp]*_Tc*Utility::pow<3>(_mag_y[_qp]) + 3*_test[_i][_qp]*_Tc*_mag_y[_qp]*Utility::pow<2>(_mag_z[_qp]) - 5*_test[_i][_qp]*_Tc*_mag_y[_qp] + 5*_test[_i][_qp]*_temperature[_qp]*_mag_y[_qp])/(5*_Tc - 5*_temperature[_qp]));
      }
      else if (_component == 2)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((3*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_x[_qp])*_mag_z[_qp] + 3*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_y[_qp])*_mag_z[_qp] + 3*_test[_i][_qp]*_Tc*Utility::pow<3>(_mag_z[_qp]) - 5*_test[_i][_qp]*_Tc*_mag_z[_qp] + 5*_test[_i][_qp]*_temperature[_qp]*_mag_z[_qp])/(5*_Tc - 5*_temperature[_qp]));
      }
      else
        return 0.0;
    }
  }
  if (!_alpha_long)
    mooseError("form = penalty needs the material property alpha_long");
   Real _temp1 = Utility::pow<2>(_mag_x[_qp])+Utility::pow<2>(_mag_y[_qp])+Utility::pow<2>(_mag_z[_qp]);
   Real _temp=_g0*_alpha_long[_qp]*_temp1*(_temp1-_fD)/(1.0+Utility::pow<2>(_alpha[_qp]));
  if (_component == 0)
  {
   return _temp*_mag_x[_qp]*_test[_i][_qp];
  }
  else if (_component == 1)
  {
   return _temp*_mag_y[_qp]*_test[_i][_qp];
  }
  else if (_component == 2)
  {
   return _temp*_mag_z[_qp]*_test[_i][_qp];
  }
  else
    return 0.0;
}

Real
MasterLongitudinalLLB::computeQpJacobian()
{
  if (_llb)
  {
    if (_temperature[_qp] < _Tc)
    {
      if (_component == 0)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((1.0/2.0)*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_me[_qp]) - 3*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_x[_qp]) - _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_y[_qp]) - _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_z[_qp]))/Utility::pow<2>(_me[_qp]));
      }
      else if (_component == 1)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((1.0/2.0)*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_me[_qp]) - _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_x[_qp]) - 3*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_y[_qp]) - _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_z[_qp]))/Utility::pow<2>(_me[_qp]));
      }
      else if (_component == 2)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((1.0/2.0)*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_me[_qp]) - _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_x[_qp]) - _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_y[_qp]) - 3*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(_mag_z[_qp]))/Utility::pow<2>(_me[_qp]));
      }
      else
        return 0.0;
    }
    else
    {
      if (_component == 0)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((9*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_x[_qp]) + 3*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_y[_qp]) + 3*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_z[_qp]) - 5*_phi[_j][_qp]*_test[_i][_qp]*_Tc + 5*_phi[_j][_qp]*_test[_i][_qp]*_temperature[_qp])/(5*_Tc - 5*_temperature[_qp]));
      }
      else if (_component == 1)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((3*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_x[_qp]) + 9*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_y[_qp]) + 3*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_z[_qp]) - 5*_phi[_j][_qp]*_test[_i][_qp]*_Tc + 5*_phi[_j][_qp]*_test[_i][_qp]*_temperature[_qp])/(5*_Tc - 5*_temperature[_qp]));
      }
      else if (_component == 2)
      {
        return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*((3*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_x[_qp]) + 3*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_y[_qp]) + 9*_phi[_j][_qp]*_test[_i][_qp]*_Tc*Utility::pow<2>(_mag_z[_qp]) - 5*_phi[_j][_qp]*_test[_i][_qp]*_Tc + 5*_phi[_j][_qp]*_test[_i][_qp]*_temperature[_qp])/(5*_Tc - 5*_temperature[_qp]));
      }
      else
        return 0.0;
    }
  }
  if (_component == 0)
  {
  return (_alpha_long[_qp]*_g0*(5.0*Utility::pow<4>(_mag_x[_qp]) + 6.0*Utility::pow<2>(_mag_x[_qp])*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) - _fD*(3.0*Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*
     _phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 1)
  {
  return (_alpha_long[_qp]*_g0*(Utility::pow<4>(_mag_x[_qp]) + 5.0*Utility::pow<4>(_mag_y[_qp]) + 6.0*Utility::pow<2>(_mag_y[_qp])*Utility::pow<2>(_mag_z[_qp]) + Utility::pow<4>(_mag_z[_qp]) + 2*Utility::pow<2>(_mag_x[_qp])*(3.0*Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) -
       _fD*(Utility::pow<2>(_mag_x[_qp]) + 3.0*Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 2)
  {
  return (_alpha_long[_qp]*_g0*(Utility::pow<4>(_mag_x[_qp]) + Utility::pow<4>(_mag_y[_qp]) + 6.0*Utility::pow<2>(_mag_y[_qp])*Utility::pow<2>(_mag_z[_qp]) + 5.0*Utility::pow<4>(_mag_z[_qp]) + 2*Utility::pow<2>(_mag_x[_qp])*(Utility::pow<2>(_mag_y[_qp]) + 3.0*Utility::pow<2>(_mag_z[_qp])) -
       _fD*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + 3.0*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
  }
  else
    return 0.0;
}

Real
MasterLongitudinalLLB::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_llb)
  {
    if (_temperature[_qp] < _Tc)
    {
      if (_component == 0)
      {
        if (jvar == _mag_y_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(-_phi[_j][_qp]*_test[_i][_qp]*_mag_x[_qp]*_mag_y[_qp]/Utility::pow<2>(_me[_qp]));
        }
        else if (jvar == _mag_z_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(-_phi[_j][_qp]*_test[_i][_qp]*_mag_x[_qp]*_mag_z[_qp]/Utility::pow<2>(_me[_qp]));
        }
        else
          return 0.0;
      }
      else if (_component == 1)
      {
        if (jvar == _mag_x_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(-_phi[_j][_qp]*_test[_i][_qp]*_mag_x[_qp]*_mag_y[_qp]/Utility::pow<2>(_me[_qp]));
        }
        else if (jvar == _mag_z_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(-_phi[_j][_qp]*_test[_i][_qp]*_mag_y[_qp]*_mag_z[_qp]/Utility::pow<2>(_me[_qp]));
        }
        else
          return 0.0;
      }
      else if (_component == 2)
      {
        if (jvar == _mag_x_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(-_phi[_j][_qp]*_test[_i][_qp]*_mag_x[_qp]*_mag_z[_qp]/Utility::pow<2>(_me[_qp]));
        }
        else if (jvar == _mag_y_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(-_phi[_j][_qp]*_test[_i][_qp]*_mag_y[_qp]*_mag_z[_qp]/Utility::pow<2>(_me[_qp]));
        }
        else
          return 0.0;
      }
      else
        return 0.0;
    }
    else
    {
      if (_component == 0)
      {
        if (jvar == _mag_y_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(6*_phi[_j][_qp]*_test[_i][_qp]*_Tc*_mag_x[_qp]*_mag_y[_qp]/(5*_Tc - 5*_temperature[_qp]));
        }
        else if (jvar == _mag_z_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(6*_phi[_j][_qp]*_test[_i][_qp]*_Tc*_mag_x[_qp]*_mag_z[_qp]/(5*_Tc - 5*_temperature[_qp]));
        }
        else
          return 0.0;
      }
      else if (_component == 1)
      {
        if (jvar == _mag_x_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(6*_phi[_j][_qp]*_test[_i][_qp]*_Tc*_mag_x[_qp]*_mag_y[_qp]/(5*_Tc - 5*_temperature[_qp]));
        }
        else if (jvar == _mag_z_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(6*_phi[_j][_qp]*_test[_i][_qp]*_Tc*_mag_y[_qp]*_mag_z[_qp]/(5*_Tc - 5*_temperature[_qp]));
        }
        else
          return 0.0;
      }
      else if (_component == 2)
      {
        if (jvar == _mag_x_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(6*_phi[_j][_qp]*_test[_i][_qp]*_Tc*_mag_x[_qp]*_mag_z[_qp]/(5*_Tc - 5*_temperature[_qp]));
        }
        else if (jvar == _mag_y_var)
        {
          return -_alpha_par[_qp]*_g0/(_Hscale*_chi_par[_qp]*(Utility::pow<2>(_alpha[_qp]) + 1))*(6*_phi[_j][_qp]*_test[_i][_qp]*_Tc*_mag_y[_qp]*_mag_z[_qp]/(5*_Tc - 5*_temperature[_qp]));
        }
        else
          return 0.0;
      }
      else
        return 0.0;
    }
  }
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
    return (2.0*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_y[_qp]*(-_fD + 2*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
    return (2.0*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_z[_qp]*(-_fD + 2*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
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
    return (2.0*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_y[_qp]*(-_fD + 2*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
    return (2.0*_alpha_long[_qp]*_g0*_mag_y[_qp]*_mag_z[_qp]*(-_fD + 2*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
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
    return (2.0*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_z[_qp]*(-_fD + 2*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
    return (2.0*_alpha_long[_qp]*_g0*_mag_y[_qp]*_mag_z[_qp]*(-_fD + 2*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
