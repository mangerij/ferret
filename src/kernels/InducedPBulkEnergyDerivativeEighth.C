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

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "InducedPBulkEnergyDerivativeEighth.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", InducedPBulkEnergyDerivativeEighth);

template<>
InputParameters validParams<InducedPBulkEnergyDerivativeEighth>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2.0 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("induced_polar_x", "The x component of the induced polarization");
  params.addRequiredCoupledVar("induced_polar_y", "The y component of the induced polarization");
  params.addCoupledVar("induced_polar_z", 0.0, "The z component of the induced polarization");
  return params;
}

InducedPBulkEnergyDerivativeEighth::InducedPBulkEnergyDerivativeEighth(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _induced_polar_x_var(coupled("induced_polar_x")),
   _induced_polar_y_var(coupled("induced_polar_y")),
   _induced_polar_z_var(coupled("induced_polar_z")),
   _induced_polar_x(coupledValue("induced_polar_x")),
   _induced_polar_y(coupledValue("induced_polar_y")),
   _induced_polar_z(coupledValue("induced_polar_z")),
   _alpha1(getMaterialProperty<Real>("alpha1")),
   _alpha11(getMaterialProperty<Real>("alpha11")),
   _alpha12(getMaterialProperty<Real>("alpha12")),
   _alpha111(getMaterialProperty<Real>("alpha111")),
   _alpha112(getMaterialProperty<Real>("alpha112")),
   _alpha123(getMaterialProperty<Real>("alpha123")),
   _alpha1111(getMaterialProperty<Real>("alpha1111")),
   _alpha1112(getMaterialProperty<Real>("alpha1112")),
   _alpha1122(getMaterialProperty<Real>("alpha1122")),
   _alpha1123(getMaterialProperty<Real>("alpha1123"))
{
}

Real
InducedPBulkEnergyDerivativeEighth::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (0. + 6.0*_alpha111[_qp]*Utility::pow<5>(_induced_polar_x[_qp]) + 2.0*_alpha12[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp]) + 2.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<4>(_induced_polar_y[_qp]) + 2.0*_alpha12[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + 
   2.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<4>(_induced_polar_z[_qp]) + 30*_alpha111[_qp]*Utility::pow<4>(_induced_polar_x[_qp])*_polar_x[_qp] + 2.0*_alpha12[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_y[_qp])*_polar_x[_qp] + 2.0*_alpha12[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp] + 
   2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_z[_qp])*_polar_x[_qp] + 60*_alpha111[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 60*_alpha111[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<3>(_polar_x[_qp]) + 
   4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<3>(_polar_x[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<3>(_polar_x[_qp]) + 30*_alpha111[_qp]*_induced_polar_x[_qp]*Utility::pow<4>(_polar_x[_qp]) + 6.0*_alpha111[_qp]*Utility::pow<5>(_polar_x[_qp]) + 2.0*_alpha1[_qp]*(_induced_polar_x[_qp] + _polar_x[_qp]) + 4.0*_alpha11[_qp]*Utility::pow<3>(_induced_polar_x[_qp] + _polar_x[_qp]) + 4.0*_alpha12[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_y[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_polar_y[_qp] + 
   8*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_polar_y[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp] + 4.0*_alpha12[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + 24.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_polar_x[_qp]*_polar_y[_qp] + 4.0*_alpha123[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 
   8*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + 2.0*_alpha12[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2.0*_alpha12[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 
   12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 
   8*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + 8*_alpha112[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]) + 2.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<4>(_polar_y[_qp]) + 2.0*_alpha112[_qp]*_polar_x[_qp]*Utility::pow<4>(_polar_y[_qp]) + 4.0*_alpha12[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*_polar_z[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 
   8*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_z[_qp] + 4.0*_alpha12[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 24.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 4.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_x[_qp]*_polar_z[_qp] + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 8*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + 
   8*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + 8*_alpha123[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 2.0*_alpha12[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
   12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha12[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
   12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha123[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
   8*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 8*_alpha112[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]) + 2.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<4>(_polar_z[_qp]) + 2.0*_alpha112[_qp]*_polar_x[_qp]*Utility::pow<4>(_polar_z[_qp]));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (0. + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_x[_qp])*_induced_polar_y[_qp] + 4.0*_alpha11[_qp]*Utility::pow<3>(_induced_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<3>(_induced_polar_y[_qp]) + 6.0*_alpha111[_qp]*Utility::pow<5>(_induced_polar_y[_qp]) + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + 2.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<4>(_induced_polar_z[_qp]) + 
   8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_polar_x[_qp] + 8*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_polar_x[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 
   8*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_polar_x[_qp]) + 2.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<4>(_polar_x[_qp]) + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_x[_qp])*_polar_y[_qp] + 12.0*_alpha11[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_y[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp])*_polar_y[_qp] + 30*_alpha111[_qp]*Utility::pow<4>(_induced_polar_y[_qp])*_polar_y[_qp] + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp] + 
   12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_z[_qp])*_polar_y[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp]*_polar_y[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 
   12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 8*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_polar_x[_qp])*_polar_y[_qp] + 12.0*_alpha11[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*Utility::pow<2>(_polar_y[_qp]) + 
   60*_alpha111[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 4.0*_alpha11[_qp]*Utility::pow<3>(_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp]) + 
   60*_alpha111[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<3>(_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<3>(_polar_y[_qp]) + 8*_alpha112[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp]) + 30*_alpha111[_qp]*_induced_polar_y[_qp]*Utility::pow<4>(_polar_y[_qp]) + 6.0*_alpha111[_qp]*Utility::pow<5>(_polar_y[_qp]) + 2.0*_alpha1[_qp]*(_induced_polar_y[_qp] + _polar_y[_qp]) + 
   4.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_z[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 8*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_z[_qp] + 8*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 4.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + 
   24.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_y[_qp]*_polar_z[_qp] + 8*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + 24.0*_alpha112[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 8*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + 
   2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
   2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
   12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 8*_alpha112[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 8*_alpha112[_qp]*_induced_polar_z[_qp]*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]) + 2.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<4>(_polar_z[_qp]) + 2.0*_alpha112[_qp]*_polar_y[_qp]*Utility::pow<4>(_polar_z[_qp]) + 
   2.0*_alpha12[_qp]*(_induced_polar_y[_qp] + _polar_y[_qp])*(Utility::pow<2>(_induced_polar_x[_qp] + _polar_x[_qp]) + Utility::pow<2>(_induced_polar_z[_qp] + _polar_z[_qp])));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (0. + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_x[_qp])*_induced_polar_z[_qp] + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_y[_qp])*_induced_polar_z[_qp] + 4.0*_alpha11[_qp]*Utility::pow<3>(_induced_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<3>(_induced_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<3>(_induced_polar_z[_qp]) + 6.0*_alpha111[_qp]*Utility::pow<5>(_induced_polar_z[_qp]) + 
   8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*_polar_x[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_x[_qp] + 8*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 
   8*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_x[_qp]) + 2.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<4>(_polar_x[_qp]) + 4.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_y[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_y[_qp] + 8*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_y[_qp] + 8*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 4.0*_alpha123[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 
   2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 8*_alpha112[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_y[_qp]) + 
   2.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<4>(_polar_y[_qp]) + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_x[_qp])*_polar_z[_qp] + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp])*_polar_z[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_induced_polar_y[_qp])*_polar_z[_qp] + 12.0*_alpha11[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_z[_qp])*_polar_z[_qp] + 
   12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp])*_polar_z[_qp] + 30*_alpha111[_qp]*Utility::pow<4>(_induced_polar_z[_qp])*_polar_z[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 
   2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 8*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_polar_x[_qp])*_polar_z[_qp] + 4.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_polar_y[_qp]*_polar_z[_qp] + 8*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 
   24.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp]*_polar_z[_qp] + 8*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 
   4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 2.0*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 8*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + 2.0*_alpha112[_qp]*Utility::pow<4>(_polar_y[_qp])*_polar_z[_qp] + 12.0*_alpha11[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
   12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*Utility::pow<2>(_polar_z[_qp]) + 60*_alpha111[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 24.0*_alpha112[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
   12.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha11[_qp]*Utility::pow<3>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<3>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]) + 60*_alpha111[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<3>(_polar_z[_qp]) + 8*_alpha112[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]) + 
   4.0*_alpha112[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<3>(_polar_z[_qp]) + 8*_alpha112[_qp]*_induced_polar_y[_qp]*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]) + 30*_alpha111[_qp]*_induced_polar_z[_qp]*Utility::pow<4>(_polar_z[_qp]) + 6.0*_alpha111[_qp]*Utility::pow<5>(_polar_z[_qp]) + 2.0*_alpha1[_qp]*(_induced_polar_z[_qp] + _polar_z[_qp]) + 
   2.0*_alpha12[_qp]*(Utility::pow<2>(_induced_polar_x[_qp] + _polar_x[_qp]) + Utility::pow<2>(_induced_polar_y[_qp] + _polar_y[_qp]))*(_induced_polar_z[_qp] + _polar_z[_qp]));
  }
  else
    return 0.0;
}

Real
InducedPBulkEnergyDerivativeEighth::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2.0*(_alpha1[_qp] + 15.0*_alpha111[_qp]*Utility::pow<4>(_induced_polar_x[_qp]) + _alpha12[_qp]*Utility::pow<2>(_induced_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp]) + _alpha112[_qp]*Utility::pow<4>(_induced_polar_y[_qp]) + _alpha12[_qp]*Utility::pow<2>(_induced_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + _alpha112[_qp]*Utility::pow<4>(_induced_polar_z[_qp]) + 
     60*_alpha111[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp] + 90.0*_alpha111[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 
     60*_alpha111[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + 15.0*_alpha111[_qp]*Utility::pow<4>(_polar_x[_qp]) + 6.0*_alpha11[_qp]*Utility::pow<2>(_induced_polar_x[_qp] + _polar_x[_qp]) + 2.0*_alpha12[_qp]*_induced_polar_y[_qp]*_polar_y[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_polar_y[_qp] + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_polar_y[_qp] + 2.0*_alpha123[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp] + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + 
     12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _alpha12[_qp]*Utility::pow<2>(_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 
     4.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + _alpha112[_qp]*Utility::pow<4>(_polar_y[_qp]) + 2.0*_alpha12[_qp]*_induced_polar_z[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_z[_qp] + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 
     4.0*_alpha123[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + 2.0*_alpha123[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _alpha12[_qp]*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
     6.0*_alpha112[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + _alpha123[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + _alpha112[_qp]*Utility::pow<4>(_polar_z[_qp])));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2.0*(_alpha1[_qp] + _alpha112[_qp]*Utility::pow<4>(_induced_polar_x[_qp]) + 6.0*_alpha11[_qp]*Utility::pow<2>(_induced_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp]) + 15.0*_alpha111[_qp]*Utility::pow<4>(_induced_polar_y[_qp]) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + _alpha112[_qp]*Utility::pow<4>(_induced_polar_z[_qp]) + 
     4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp] + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp] + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 4.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + 
     _alpha112[_qp]*Utility::pow<4>(_polar_x[_qp]) + 12.0*_alpha11[_qp]*_induced_polar_y[_qp]*_polar_y[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_polar_y[_qp] + 60*_alpha111[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_polar_y[_qp] + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp] + 24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 6.0*_alpha11[_qp]*Utility::pow<2>(_polar_y[_qp]) + 
     6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 90.0*_alpha111[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 60*_alpha111[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + 
     15.0*_alpha111[_qp]*Utility::pow<4>(_polar_y[_qp]) + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_z[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2.0*_alpha123[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 24.0*_alpha112[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 
     _alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + _alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
     6.0*_alpha112[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 4.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + _alpha112[_qp]*Utility::pow<4>(_polar_z[_qp]) + _alpha12[_qp]*(Utility::pow<2>(_induced_polar_x[_qp]) + Utility::pow<2>(_induced_polar_z[_qp]) + 2.0*_induced_polar_x[_qp]*_polar_x[_qp] + Utility::pow<2>(_polar_x[_qp]) + 2.0*_induced_polar_z[_qp]*_polar_z[_qp] + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2.0*(_alpha1[_qp] + _alpha112[_qp]*Utility::pow<4>(_induced_polar_x[_qp]) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_y[_qp]) + _alpha112[_qp]*Utility::pow<4>(_induced_polar_y[_qp]) + 6.0*_alpha11[_qp]*Utility::pow<2>(_induced_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_induced_polar_z[_qp]) + 15.0*_alpha111[_qp]*Utility::pow<4>(_induced_polar_z[_qp]) + 
     4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_x[_qp])*_polar_x[_qp] + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_polar_x[_qp] + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_x[_qp] + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 4.0*_alpha112[_qp]*_induced_polar_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + 
     _alpha112[_qp]*Utility::pow<4>(_polar_x[_qp]) + 2.0*_alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_y[_qp]*_polar_y[_qp] + 4.0*_alpha112[_qp]*Utility::pow<3>(_induced_polar_y[_qp])*_polar_y[_qp] + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*_polar_y[_qp] + 4.0*_alpha123[_qp]*_induced_polar_x[_qp]*_induced_polar_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2.0*_alpha123[_qp]*_induced_polar_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 
     6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2.0*_alpha123[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + _alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 4.0*_alpha112[_qp]*_induced_polar_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + _alpha112[_qp]*Utility::pow<4>(_polar_y[_qp]) + 
     _alpha12[_qp]*(Utility::pow<2>(_induced_polar_x[_qp]) + Utility::pow<2>(_induced_polar_y[_qp]) + 2.0*_induced_polar_x[_qp]*_polar_x[_qp] + Utility::pow<2>(_polar_x[_qp]) + 2.0*_induced_polar_y[_qp]*_polar_y[_qp] + Utility::pow<2>(_polar_y[_qp])) + 12.0*_alpha11[_qp]*_induced_polar_z[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*_induced_polar_z[_qp]*_polar_z[_qp] + 60*_alpha111[_qp]*Utility::pow<3>(_induced_polar_z[_qp])*_polar_z[_qp] + 
     24.0*_alpha112[_qp]*_induced_polar_x[_qp]*_induced_polar_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 24.0*_alpha112[_qp]*_induced_polar_y[_qp]*_induced_polar_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + 12.0*_alpha112[_qp]*_induced_polar_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 6.0*_alpha11[_qp]*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_induced_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
     90.0*_alpha111[_qp]*Utility::pow<2>(_induced_polar_z[_qp])*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 12.0*_alpha112[_qp]*_induced_polar_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 6.0*_alpha112[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 60*_alpha111[_qp]*_induced_polar_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 15.0*_alpha111[_qp]*Utility::pow<4>(_polar_z[_qp])));
  }
  else
    return 0.0;
}

Real
InducedPBulkEnergyDerivativeEighth::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _induced_polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4.0*(_induced_polar_x[_qp] + _polar_x[_qp])*(_induced_polar_y[_qp] + _polar_y[_qp])*(_alpha12[_qp] + 2.0*_alpha112[_qp]*(Utility::pow<2>(_induced_polar_x[_qp]) + Utility::pow<2>(_induced_polar_y[_qp]) + 2.0*_induced_polar_x[_qp]*_polar_x[_qp] + Utility::pow<2>(_polar_x[_qp]) + 2.0*_induced_polar_y[_qp]*_polar_y[_qp] + Utility::pow<2>(_polar_y[_qp])) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_z[_qp] + _polar_z[_qp])));
    }
    else if (jvar == _induced_polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4.0*(_induced_polar_x[_qp] + _polar_x[_qp])*(_induced_polar_z[_qp] + _polar_z[_qp])*(_alpha12[_qp] + _alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp] + _polar_y[_qp]) + 2.0*_alpha112[_qp]*(Utility::pow<2>(_induced_polar_x[_qp]) + Utility::pow<2>(_induced_polar_z[_qp]) + 2.0*_induced_polar_x[_qp]*_polar_x[_qp] + Utility::pow<2>(_polar_x[_qp]) + 2.0*_induced_polar_z[_qp]*_polar_z[_qp] + Utility::pow<2>(_polar_z[_qp]))));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _induced_polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4.0*(_induced_polar_x[_qp] + _polar_x[_qp])*(_induced_polar_y[_qp] + _polar_y[_qp])*(_alpha12[_qp] + 2.0*_alpha112[_qp]*(Utility::pow<2>(_induced_polar_x[_qp]) + Utility::pow<2>(_induced_polar_y[_qp]) + 2.0*_induced_polar_x[_qp]*_polar_x[_qp] + Utility::pow<2>(_polar_x[_qp]) + 2.0*_induced_polar_y[_qp]*_polar_y[_qp] + Utility::pow<2>(_polar_y[_qp])) + _alpha123[_qp]*Utility::pow<2>(_induced_polar_z[_qp] + _polar_z[_qp])));
    }
    else if (jvar == _induced_polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4.0*(_induced_polar_y[_qp] + _polar_y[_qp])*(_induced_polar_z[_qp] + _polar_z[_qp])*(_alpha12[_qp] + _alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp] + _polar_x[_qp]) + 2.0*_alpha112[_qp]*(Utility::pow<2>(_induced_polar_y[_qp]) + Utility::pow<2>(_induced_polar_z[_qp]) + 2.0*_induced_polar_y[_qp]*_polar_y[_qp] + Utility::pow<2>(_polar_y[_qp]) + 2.0*_induced_polar_z[_qp]*_polar_z[_qp] + Utility::pow<2>(_polar_z[_qp]))));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _induced_polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4.0*(_induced_polar_x[_qp] + _polar_x[_qp])*(_induced_polar_z[_qp] + _polar_z[_qp])*(_alpha12[_qp] + _alpha123[_qp]*Utility::pow<2>(_induced_polar_y[_qp] + _polar_y[_qp]) + 2.0*_alpha112[_qp]*(Utility::pow<2>(_induced_polar_x[_qp]) + Utility::pow<2>(_induced_polar_z[_qp]) + 2.0*_induced_polar_x[_qp]*_polar_x[_qp] + Utility::pow<2>(_polar_x[_qp]) + 2.0*_induced_polar_z[_qp]*_polar_z[_qp] + Utility::pow<2>(_polar_z[_qp]))));
    }
    else if (jvar == _induced_polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4.0*(_induced_polar_y[_qp] + _polar_y[_qp])*(_induced_polar_z[_qp] + _polar_z[_qp])*(_alpha12[_qp] + _alpha123[_qp]*Utility::pow<2>(_induced_polar_x[_qp] + _polar_x[_qp]) + 2.0*_alpha112[_qp]*(Utility::pow<2>(_induced_polar_y[_qp]) + Utility::pow<2>(_induced_polar_z[_qp]) + 2.0*_induced_polar_y[_qp]*_polar_y[_qp] + Utility::pow<2>(_polar_y[_qp]) + 2.0*_induced_polar_z[_qp]*_polar_z[_qp] + Utility::pow<2>(_polar_z[_qp]))));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
