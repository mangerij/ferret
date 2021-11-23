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

#include "RotoBulkEnergyDerivativeEighthAlt.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotoBulkEnergyDerivativeEighthAlt);

InputParameters RotoBulkEnergyDerivativeEighthAlt::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotoBulkEnergyDerivativeEighthAlt::RotoBulkEnergyDerivativeEighthAlt(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _beta1(getMaterialProperty<Real>("beta1")),
   _beta11(getMaterialProperty<Real>("beta11")),
   _beta12(getMaterialProperty<Real>("beta12")),
   _beta111(getMaterialProperty<Real>("beta111")),
   _beta112(getMaterialProperty<Real>("beta112")),
   _beta123(getMaterialProperty<Real>("beta123")),
   _beta1111(getMaterialProperty<Real>("beta1111")),
   _beta1112(getMaterialProperty<Real>("beta1112")),
   _beta1122(getMaterialProperty<Real>("beta1122")),
   _beta1123(getMaterialProperty<Real>("beta1123"))
{
}

Real
RotoBulkEnergyDerivativeEighthAlt::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_beta1[_qp]*_antiferrodis_A_x[_qp] + 4*_beta11[_qp]*Utility::pow<3>(_antiferrodis_A_x[_qp]) + 6*_beta111[_qp]*Utility::pow<5>(_antiferrodis_A_x[_qp]) + 8*_beta1111[_qp]*Utility::pow<7>(_antiferrodis_A_x[_qp]) + 
   2*_beta123[_qp]*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12[_qp]*(2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1122[_qp]*(4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112[_qp]*(2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      4*Utility::pow<3>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112[_qp]*(2*_antiferrodis_A_x[_qp]*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      6*Utility::pow<5>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*_beta1[_qp]*_antiferrodis_A_y[_qp] + 4*_beta11[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 6*_beta111[_qp]*Utility::pow<5>(_antiferrodis_A_y[_qp]) + 8*_beta1111[_qp]*Utility::pow<7>(_antiferrodis_A_y[_qp]) + 
   2*_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12[_qp]*(2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 4*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1122[_qp]*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      4*Utility::pow<3>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 2*_antiferrodis_A_y[_qp]*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      6*Utility::pow<5>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_beta1[_qp]*_antiferrodis_A_z[_qp] + 2*_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 4*_beta11[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   6*_beta111[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp]) + 8*_beta1111[_qp]*Utility::pow<7>(_antiferrodis_A_z[_qp]) + 
   _beta12[_qp]*(2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]) + 
   _beta1123[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      4*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1122[_qp]*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 4*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      4*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 2*Utility::pow<6>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      6*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<5>(_antiferrodis_A_z[_qp])));
  }
  else
    return 0.0;
}

Real
RotoBulkEnergyDerivativeEighthAlt::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1[_qp] + 12*_beta11[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 30*_beta111[_qp]*Utility::pow<4>(_antiferrodis_A_x[_qp]) + 56*_beta1111[_qp]*Utility::pow<6>(_antiferrodis_A_x[_qp]) + 
   2*_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12[_qp]*(2*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1122[_qp]*(12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      12*Utility::pow<2>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      30*Utility::pow<4>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1[_qp] + 12*_beta11[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 30*_beta111[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 56*_beta1111[_qp]*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 
   2*_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + _beta12[_qp]*(2*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta1122[_qp]*(12*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_antiferrodis_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 
      12*Utility::pow<2>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
      30*Utility::pow<4>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_beta1[_qp] + 2*_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp]) + _beta12[_qp]*(2*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])) + 
   12*_beta11[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 30*_beta111[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp]) + 56*_beta1111[_qp]*Utility::pow<6>(_antiferrodis_A_z[_qp]) + 
   _beta1123[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 
      12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1122[_qp]*(12*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 12*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta112[_qp]*(2*Utility::pow<4>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<4>(_antiferrodis_A_y[_qp]) + 
      12*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<2>(_antiferrodis_A_z[_qp])) + 
   _beta1112[_qp]*(2*Utility::pow<6>(_antiferrodis_A_x[_qp]) + 2*Utility::pow<6>(_antiferrodis_A_y[_qp]) + 
      30*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<4>(_antiferrodis_A_z[_qp])));
  }
  else
    return 0.0;
}

Real
RotoBulkEnergyDerivativeEighthAlt::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])) + 
   _beta1112[_qp]*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_y[_qp])) + 4*_beta123[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp])));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123[_qp]*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 4*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])) + 
   _beta1112[_qp]*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_y[_qp])) + 4*_beta123[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_antiferrodis_A_z[_qp]) + 
      4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<4>(_antiferrodis_A_z[_qp])));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 8*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_y[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123[_qp]*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(8*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 4*_antiferrodis_A_x[_qp]*Utility::pow<4>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_beta12[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 4*_beta123[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 16*_beta1122[_qp]*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_antiferrodis_A_z[_qp]) + 
   _beta112[_qp]*(8*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 8*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + 
   _beta1123[_qp]*(4*Utility::pow<4>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp] + 8*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 
      8*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<3>(_antiferrodis_A_z[_qp])) + _beta1112[_qp]*(12*Utility::pow<5>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] + 12*_antiferrodis_A_y[_qp]*Utility::pow<5>(_antiferrodis_A_z[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
