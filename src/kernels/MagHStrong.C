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

#include "MagHStrong.h"

class MagHStrong;

registerMooseObject("FerretApp", MagHStrong);

InputParameters MagHStrong::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for bound magnetic charge (div M)");
  params.addRequiredCoupledVar("azimuthal_ph", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_th", "The polar component of the constrained magnetic vector");
  return params;
}

MagHStrong::MagHStrong(const InputParameters & parameters)
  :Kernel(parameters),
   _azimuthal_ph_var(coupled("azimuthal_ph")),
   _polar_th_var(coupled("polar_th")),
   _azimuthal_ph(coupledValue("azimuthal_ph")),
   _polar_th(coupledValue("polar_th")),
   _mu0(getMaterialProperty<Real>("mu0")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
MagHStrong::computeQpResidual()
{
  return -_mu0[_qp]*(_Ms[_qp]*(_grad_test[_i][_qp](2)*std::cos(_polar_th[_qp])+(_grad_test[_i][_qp](0)*std::cos(_azimuthal_ph[_qp])+_grad_test[_i][_qp](1)*std::sin(_azimuthal_ph[_qp]))*std::sin(_polar_th[_qp])));
}
Real
MagHStrong::computeQpJacobian()
{
  return 0.0;
}

Real
MagHStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _polar_th_var)
  {
    return -_mu0[_qp]*_Ms[_qp]*_phi[_j][_qp]*(_grad_test[_i][_qp](0)*std::cos(_azimuthal_ph[_qp])*std::cos(_polar_th[_qp])+_grad_test[_i][_qp](1)*std::sin(_azimuthal_ph[_qp])*std::cos(_polar_th[_qp])-_grad_test[_i][_qp](2)*std::sin(_polar_th[_qp]));
  }
  else if (jvar == _azimuthal_ph_var)
  {
    return -_mu0[_qp]*_Ms[_qp]*_phi[_j][_qp]*((_grad_test[_i][_qp](1)*std::cos(_azimuthal_ph[_qp]))-_grad_test[_i][_qp](0)*std::sin(_azimuthal_ph[_qp]))*std::sin(_polar_th[_qp]);
  }
  else
    return 0.0;
}
