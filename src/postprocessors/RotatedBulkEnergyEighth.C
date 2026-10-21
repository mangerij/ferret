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

#include "RotatedBulkEnergyEighth.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotatedBulkEnergyEighth);

InputParameters RotatedBulkEnergyEighth::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral whose integrand is the eighth order expansion of the polarization, with the crystal axes turned by the material property twist about z (the energy of RotatedBulkEnergyDerivativeEighth).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

RotatedBulkEnergyEighth::RotatedBulkEnergyEighth(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _alpha1(getMaterialProperty<Real>("alpha1")),
   _alpha11(getMaterialProperty<Real>("alpha11")),
   _alpha12(getMaterialProperty<Real>("alpha12")),
   _alpha111(getMaterialProperty<Real>("alpha111")),
   _alpha112(getMaterialProperty<Real>("alpha112")),
   _alpha123(getMaterialProperty<Real>("alpha123")),
   _alpha1111(getMaterialProperty<Real>("alpha1111")),
   _alpha1112(getMaterialProperty<Real>("alpha1112")),
   _alpha1122(getMaterialProperty<Real>("alpha1122")),
   _alpha1123(getMaterialProperty<Real>("alpha1123")),
   _twist(getMaterialProperty<Real>("twist")),
   _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
RotatedBulkEnergyEighth::computeQpIntegral()
{
  // polarization in the crystal frame: P_c = R(twist)^T P
  const Real px = _polar_x[_qp]*std::cos(_twist[_qp]) + _polar_y[_qp]*std::sin(_twist[_qp]);
  const Real py = -_polar_x[_qp]*std::sin(_twist[_qp]) + _polar_y[_qp]*std::cos(_twist[_qp]);
  const Real pz = _polar_z[_qp];

  return _energy_scale*(_alpha1[_qp]*(Utility::pow<2>(px) + Utility::pow<2>(py) + Utility::pow<2>(pz)) +
         _alpha11[_qp]*(Utility::pow<4>(px) + Utility::pow<4>(py) + Utility::pow<4>(pz)) +
         _alpha12[_qp]*(Utility::pow<2>(px)*Utility::pow<2>(py) + Utility::pow<2>(px)*Utility::pow<2>(pz) + Utility::pow<2>(py)*Utility::pow<2>(pz)) +
         _alpha111[_qp]*(Utility::pow<6>(px) + Utility::pow<6>(py) + Utility::pow<6>(pz)) +
         _alpha112[_qp]*(Utility::pow<4>(px)*(Utility::pow<2>(py) + Utility::pow<2>(pz)) + Utility::pow<4>(py)*(Utility::pow<2>(px) + Utility::pow<2>(pz)) + Utility::pow<4>(pz)*(Utility::pow<2>(px) + Utility::pow<2>(py))) +
         _alpha123[_qp]*Utility::pow<2>(px)*Utility::pow<2>(py)*Utility::pow<2>(pz) +
         _alpha1111[_qp]*(Utility::pow<8>(px) + Utility::pow<8>(py) + Utility::pow<8>(pz)) +
         _alpha1112[_qp]*(Utility::pow<6>(px)*(Utility::pow<2>(py) + Utility::pow<2>(pz)) + Utility::pow<6>(py)*(Utility::pow<2>(px) + Utility::pow<2>(pz)) + Utility::pow<6>(pz)*(Utility::pow<2>(px) + Utility::pow<2>(py))) +
         _alpha1122[_qp]*(Utility::pow<4>(px)*Utility::pow<4>(py) + Utility::pow<4>(px)*Utility::pow<4>(pz) + Utility::pow<4>(py)*Utility::pow<4>(pz)) +
         _alpha1123[_qp]*(Utility::pow<4>(px)*Utility::pow<2>(py)*Utility::pow<2>(pz) + Utility::pow<2>(px)*Utility::pow<4>(py)*Utility::pow<2>(pz) + Utility::pow<2>(px)*Utility::pow<2>(py)*Utility::pow<4>(pz)));
}
