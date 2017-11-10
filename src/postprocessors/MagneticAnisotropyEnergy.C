/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "MagneticAnisotropyEnergy.h"

template<>
InputParameters validParams<MagneticAnisotropyEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("antiferromag_L_x", "The x component of the antiferromagnetic vector");
  params.addCoupledVar("antiferromag_L_y", 0.0, "The y component of the antiferromagnetic vector");
  params.addCoupledVar("antiferromag_L_z", 0.0, "The z component of the antiferromagnetic vector");
  params.addRequiredParam<Real>("Ku", "The constant of anisotropy");
  params.addRequiredParam<Real>("nx", "x direction of the anisotropy");
  params.addRequiredParam<Real>("ny", "y direction of the anisotropy");
  params.addRequiredParam<Real>("nz", "z direction of the anisotropy");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

MagneticAnisotropyEnergy::MagneticAnisotropyEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _antiferromag_L_x(coupledValue("antiferromag_L_x")),
  _antiferromag_L_y(coupledValue("antiferromag_L_y")),
  _antiferromag_L_z(coupledValue("antiferromag_L_z")),
  _Ku(getParam<Real>("Ku")),
  _nx(getParam<Real>("nx")),
  _ny(getParam<Real>("ny")),
  _nz(getParam<Real>("nz")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
MagneticAnisotropyEnergy::computeQpIntegral()
{
  return _Ku*std::pow(_antiferromag_L_x[_qp]*_nx + _antiferromag_L_y[_qp]*_ny + _antiferromag_L_z[_qp]*_nz,2);
}
