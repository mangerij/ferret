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

#include "AnisotropyUSLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AnisotropyUSLLG);

InputParameters AnisotropyUSLLG::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_th", "The polar angle of the magnetic vector");
  params.addRequiredCoupledVar("azimuthal_ph", "The azimuthal angle of the magnetic vector");
  return params;
}

AnisotropyUSLLG::AnisotropyUSLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_th_var(coupled("polar_th")),
  _azimuthal_ph_var(coupled("azimuthal_ph")),
  _polar_th(coupledValue("polar_th")),
  _azimuthal_ph(coupledValue("azimuthal_ph")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _K1(getMaterialProperty<Real>("K1")),
   //_K2(getMaterialProperty<Real>("K2")),
  _nx(getMaterialProperty<Real>("nx")),
  _ny(getMaterialProperty<Real>("ny")),
  _nz(getMaterialProperty<Real>("nz")),
  _g0(getMaterialProperty<Real>("g0")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _mu0(getMaterialProperty<Real>("mu0"))
{
}

Real
AnisotropyUSLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (2.*_g0[_qp]*_K1[_qp]*_test[_i][_qp]*(std::cos(_azimuthal_ph[_qp])*(_ny[_qp] + _alpha[_qp]*_nx[_qp]*std::cos(_polar_th[_qp]))
	    + (-_nx[_qp] + _alpha[_qp]*_ny[_qp]*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp])
	       -_alpha[_qp]*_nz[_qp]*std::sin(_polar_th[_qp]))*(_nz[_qp]*std::cos(_polar_th[_qp])
	       +(_nx[_qp]*std::cos(_azimuthal_ph[_qp]) + _ny[_qp]*std::sin(_azimuthal_ph[_qp]))*std::sin(_polar_th[_qp])))
      /((1. +Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return (-2.*_g0[_qp]*_K1[_qp]*_test[_i][_qp]*(_nx[_qp]*std::cos(_azimuthal_ph[_qp]) + _nz[_qp]*std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]) + _ny[_qp]*std::sin(_azimuthal_ph[_qp]))
	     *(std::cos(_azimuthal_ph[_qp])*(-(_alpha[_qp]*_ny[_qp]) + _nx[_qp]*std::cos(_polar_th[_qp])) + (_alpha[_qp]*_nx[_qp] + _ny[_qp]*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp]) - _nz[_qp]*std::sin(_polar_th[_qp])))
      /((1. + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp])
      ;
  }
  else
    return 0.0;
}

Real
AnisotropyUSLLG::computeQpJacobian()
{
  if (_component == 0)
    {
    return (2.*_g0[_qp]*_K1[_qp]*_test[_i][_qp]*(_azimuthal_ph[_qp]*(_nx[_qp]*std::cos(_azimuthal_ph[_qp])*std::cos(_polar_th[_qp]) + _ny[_qp]*std::cos(_polar_th[_qp])*std::sin(_azimuthal_ph[_qp]) - _nz[_qp]*std::sin(_polar_th[_qp]))*(std::cos(_azimuthal_ph[_qp])*(_ny[_qp] + _alpha[_qp]*_nx[_qp]*std::cos(_polar_th[_qp])) + (-_nx[_qp] + _alpha[_qp]*_ny[_qp]*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp]) - _alpha[_qp]*_nz[_qp]*std::sin(_polar_th[_qp])) - _alpha[_qp]*_azimuthal_ph[_qp]*Utility::pow<2>(_nz[_qp]*std::cos(_polar_th[_qp]) + (_nx[_qp]*std::cos(_azimuthal_ph[_qp]) + _ny[_qp]*std::sin(_azimuthal_ph[_qp]))*std::sin(_polar_th[_qp]))))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
   else if (_component == 1)
     {
       return (2*_g0[_qp]*_K1[_qp]*_test[_i][_qp]*(-(_azimuthal_ph[_qp]*(_nx[_qp]*std::cos(_azimuthal_ph[_qp]) + _nz[_qp]*std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]) + _ny[_qp]*std::sin(_azimuthal_ph[_qp]))*(std::cos(_azimuthal_ph[_qp])*(_alpha[_qp]*_nx[_qp] + _ny[_qp]*std::cos(_polar_th[_qp])) + (_alpha[_qp]*_ny[_qp] - _nx[_qp]*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp]))) + _azimuthal_ph[_qp]*(_ny[_qp]*std::cos(_azimuthal_ph[_qp]) - _nx[_qp]*std::sin(_azimuthal_ph[_qp]))*(std::cos(_azimuthal_ph[_qp])*(_alpha[_qp]*_ny[_qp] - _nx[_qp]*std::cos(_polar_th[_qp])) - (_alpha[_qp]*_nx[_qp] + _ny[_qp]*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp]) + _nz[_qp]*std::sin(_polar_th[_qp]))))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
     }
   else
     return 0.0;
}

  Real
AnisotropyUSLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuthal_ph_var)
    {
      return (2.*_g0[_qp]*_K1[_qp]*_azimuthal_ph[_qp]*_test[_i][_qp]*(_alpha[_qp]*_nz[_qp]*std::cos(2*_polar_th[_qp])*(_ny[_qp]*std::cos(_azimuthal_ph[_qp]) - _nx[_qp]*std::sin(_azimuthal_ph[_qp])) - _nz[_qp]*std::cos(_polar_th[_qp])*(_nx[_qp]*std::cos(_azimuthal_ph[_qp]) + _ny[_qp]*std::sin(_azimuthal_ph[_qp])) + (std::cos(2*_azimuthal_ph[_qp])*(-Utility::pow<2>(_nx[_qp]) + Utility::pow<2>(_ny[_qp]) + 2*_alpha[_qp]*_nx[_qp]*_ny[_qp]*std::cos(_polar_th[_qp])) + (-2*_nx[_qp]*_ny[_qp] + _alpha[_qp]*(-Utility::pow<2>(_nx[_qp]) + Utility::pow<2>(_ny[_qp]))*std::cos(_polar_th[_qp]))*std::sin(2*_azimuthal_ph[_qp]))*std::sin(_polar_th[_qp])))/
	((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_th_var)
    {
      return (2.*_g0[_qp]*_K1[_qp]*_test[_i][_qp]*(_nz[_qp]*_azimuthal_ph[_qp]*Utility::pow<2>(1./std::sin(_polar_th[_qp]))*(std::cos(_azimuthal_ph[_qp])*(-(_alpha[_qp]*_ny[_qp]) + _nx[_qp]*std::cos(_polar_th[_qp])) + (_alpha[_qp]*_nx[_qp] + _ny[_qp]*std::cos(_polar_th[_qp]))*std::sin(_azimuthal_ph[_qp]) - _nz[_qp]*std::sin(_polar_th[_qp])) + _azimuthal_ph[_qp]*(_nx[_qp]*std::cos(_azimuthal_ph[_qp]) + _nz[_qp]*std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]) + _ny[_qp]*std::sin(_azimuthal_ph[_qp]))*(_nz[_qp]*std::cos(_polar_th[_qp]) + (_nx[_qp]*std::cos(_azimuthal_ph[_qp]) + _ny[_qp]*std::sin(_azimuthal_ph[_qp]))*std::sin(_polar_th[_qp]))))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
      return 0.0;
    }
  else
    return 0.0;
}
