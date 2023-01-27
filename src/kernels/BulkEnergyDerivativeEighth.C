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

#include "BulkEnergyDerivativeEighth.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", BulkEnergyDerivativeEighth);

InputParameters BulkEnergyDerivativeEighth::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

BulkEnergyDerivativeEighth::BulkEnergyDerivativeEighth(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
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
   _alpha1123(getMaterialProperty<Real>("alpha1123"))
{
}

Real
BulkEnergyDerivativeEighth::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_alpha1[_qp]*_polar_x[_qp] + 4*_alpha11[_qp]*Utility::pow<3>(_polar_x[_qp]) + 6*_alpha111[_qp]*Utility::pow<5>(_polar_x[_qp]) + 8*_alpha1111[_qp]*Utility::pow<7>(_polar_x[_qp]) + 
   2*_alpha123[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + _alpha12[_qp]*(2*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha1122[_qp]*(4*Utility::pow<3>(_polar_x[_qp])*Utility::pow<4>(_polar_y[_qp]) + 4*Utility::pow<3>(_polar_x[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(4*Utility::pow<3>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 2*_polar_x[_qp]*Utility::pow<4>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      2*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha112[_qp]*(2*_polar_x[_qp]*Utility::pow<4>(_polar_y[_qp]) + 2*_polar_x[_qp]*Utility::pow<4>(_polar_z[_qp]) + 
      4*Utility::pow<3>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))) + 
   _alpha1112[_qp]*(2*_polar_x[_qp]*Utility::pow<6>(_polar_y[_qp]) + 2*_polar_x[_qp]*Utility::pow<6>(_polar_z[_qp]) + 
      6*Utility::pow<5>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*_alpha1[_qp]*_polar_y[_qp] + 4*_alpha11[_qp]*Utility::pow<3>(_polar_y[_qp]) + 6*_alpha111[_qp]*Utility::pow<5>(_polar_y[_qp]) + 8*_alpha1111[_qp]*Utility::pow<7>(_polar_y[_qp]) + 
   2*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + _alpha12[_qp]*(2*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 2*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(2*Utility::pow<4>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 4*Utility::pow<2>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      2*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha1122[_qp]*(4*Utility::pow<4>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp]) + 4*Utility::pow<3>(_polar_y[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha112[_qp]*(2*Utility::pow<4>(_polar_x[_qp])*_polar_y[_qp] + 2*_polar_y[_qp]*Utility::pow<4>(_polar_z[_qp]) + 
      4*Utility::pow<3>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))) + 
   _alpha1112[_qp]*(2*Utility::pow<6>(_polar_x[_qp])*_polar_y[_qp] + 2*_polar_y[_qp]*Utility::pow<6>(_polar_z[_qp]) + 
      6*Utility::pow<5>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_alpha1[_qp]*_polar_z[_qp] + 2*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 4*_alpha11[_qp]*Utility::pow<3>(_polar_z[_qp]) + 
   6*_alpha111[_qp]*Utility::pow<5>(_polar_z[_qp]) + 8*_alpha1111[_qp]*Utility::pow<7>(_polar_z[_qp]) + 
   _alpha12[_qp]*(2*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 2*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]) + 
   _alpha1123[_qp]*(2*Utility::pow<4>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<4>(_polar_y[_qp])*_polar_z[_qp] + 
      4*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp])) + 
   _alpha1122[_qp]*(4*Utility::pow<4>(_polar_x[_qp])*Utility::pow<3>(_polar_z[_qp]) + 4*Utility::pow<4>(_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp])) + 
   _alpha112[_qp]*(2*Utility::pow<4>(_polar_x[_qp])*_polar_z[_qp] + 2*Utility::pow<4>(_polar_y[_qp])*_polar_z[_qp] + 
      4*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<3>(_polar_z[_qp])) + 
   _alpha1112[_qp]*(2*Utility::pow<6>(_polar_x[_qp])*_polar_z[_qp] + 2*Utility::pow<6>(_polar_y[_qp])*_polar_z[_qp] + 
      6*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<5>(_polar_z[_qp])));
  }
  else
    return 0.0;
}

Real
BulkEnergyDerivativeEighth::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_alpha1[_qp] + 12*_alpha11[_qp]*Utility::pow<2>(_polar_x[_qp]) + 30*_alpha111[_qp]*Utility::pow<4>(_polar_x[_qp]) + 56*_alpha1111[_qp]*Utility::pow<6>(_polar_x[_qp]) + 
   2*_alpha123[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + _alpha12[_qp]*(2*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha1122[_qp]*(12*Utility::pow<2>(_polar_x[_qp])*Utility::pow<4>(_polar_y[_qp]) + 12*Utility::pow<2>(_polar_x[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(12*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 2*Utility::pow<4>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha112[_qp]*(2*Utility::pow<4>(_polar_y[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]) + 
      12*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))) + 
   _alpha1112[_qp]*(2*Utility::pow<6>(_polar_y[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]) + 
      30*Utility::pow<4>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_alpha1[_qp] + 12*_alpha11[_qp]*Utility::pow<2>(_polar_y[_qp]) + 30*_alpha111[_qp]*Utility::pow<4>(_polar_y[_qp]) + 56*_alpha1111[_qp]*Utility::pow<6>(_polar_y[_qp]) + 
   2*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + _alpha12[_qp]*(2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(2*Utility::pow<4>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 12*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha1122[_qp]*(12*Utility::pow<4>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12*Utility::pow<2>(_polar_y[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha112[_qp]*(2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]) + 
      12*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))) + 
   _alpha1112[_qp]*(2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]) + 
      30*Utility::pow<4>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*_alpha1[_qp] + 2*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + _alpha12[_qp]*(2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp])) + 
   12*_alpha11[_qp]*Utility::pow<2>(_polar_z[_qp]) + 30*_alpha111[_qp]*Utility::pow<4>(_polar_z[_qp]) + 56*_alpha1111[_qp]*Utility::pow<6>(_polar_z[_qp]) + 
   _alpha1123[_qp]*(2*Utility::pow<4>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<4>(_polar_y[_qp]) + 
      12*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha1122[_qp]*(12*Utility::pow<4>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + 12*Utility::pow<4>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha112[_qp]*(2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_y[_qp]) + 
      12*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha1112[_qp]*(2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_y[_qp]) + 
      30*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<4>(_polar_z[_qp])));
  }
  else
    return 0.0;
}

Real
BulkEnergyDerivativeEighth::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12[_qp]*_polar_x[_qp]*_polar_y[_qp] + 16*_alpha1122[_qp]*Utility::pow<3>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp]) + 
   _alpha112[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + 8*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp])) + 
   _alpha1112[_qp]*(12*Utility::pow<5>(_polar_x[_qp])*_polar_y[_qp] + 12*_polar_x[_qp]*Utility::pow<5>(_polar_y[_qp])) + 4*_alpha123[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
   _alpha1123[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 8*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      4*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<4>(_polar_z[_qp])));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12[_qp]*_polar_x[_qp]*_polar_z[_qp] + 4*_alpha123[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 16*_alpha1122[_qp]*Utility::pow<3>(_polar_x[_qp])*Utility::pow<3>(_polar_z[_qp]) + 
   _alpha112[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + 8*_polar_x[_qp]*Utility::pow<3>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 4*_polar_x[_qp]*Utility::pow<4>(_polar_y[_qp])*_polar_z[_qp] + 
      8*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp])) + _alpha1112[_qp]*(12*Utility::pow<5>(_polar_x[_qp])*_polar_z[_qp] + 12*_polar_x[_qp]*Utility::pow<5>(_polar_z[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12[_qp]*_polar_x[_qp]*_polar_y[_qp] + 16*_alpha1122[_qp]*Utility::pow<3>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp]) + 
   _alpha112[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + 8*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp])) + 
   _alpha1112[_qp]*(12*Utility::pow<5>(_polar_x[_qp])*_polar_y[_qp] + 12*_polar_x[_qp]*Utility::pow<5>(_polar_y[_qp])) + 4*_alpha123[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 
   _alpha1123[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]) + 8*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + 
      4*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<4>(_polar_z[_qp])));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12[_qp]*_polar_y[_qp]*_polar_z[_qp] + 4*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + 16*_alpha1122[_qp]*Utility::pow<3>(_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]) + 
   _alpha112[_qp]*(8*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + 8*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(4*Utility::pow<4>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + 8*Utility::pow<2>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + 
      8*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp])) + _alpha1112[_qp]*(12*Utility::pow<5>(_polar_y[_qp])*_polar_z[_qp] + 12*_polar_y[_qp]*Utility::pow<5>(_polar_z[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12[_qp]*_polar_x[_qp]*_polar_z[_qp] + 4*_alpha123[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 16*_alpha1122[_qp]*Utility::pow<3>(_polar_x[_qp])*Utility::pow<3>(_polar_z[_qp]) + 
   _alpha112[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + 8*_polar_x[_qp]*Utility::pow<3>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(8*Utility::pow<3>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + 4*_polar_x[_qp]*Utility::pow<4>(_polar_y[_qp])*_polar_z[_qp] + 
      8*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp])) + _alpha1112[_qp]*(12*Utility::pow<5>(_polar_x[_qp])*_polar_z[_qp] + 12*_polar_x[_qp]*Utility::pow<5>(_polar_z[_qp])));
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_alpha12[_qp]*_polar_y[_qp]*_polar_z[_qp] + 4*_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + 16*_alpha1122[_qp]*Utility::pow<3>(_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]) + 
   _alpha112[_qp]*(8*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + 8*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp])) + 
   _alpha1123[_qp]*(4*Utility::pow<4>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + 8*Utility::pow<2>(_polar_x[_qp])*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + 
      8*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp])) + _alpha1112[_qp]*(12*Utility::pow<5>(_polar_y[_qp])*_polar_z[_qp] + 12*_polar_y[_qp]*Utility::pow<5>(_polar_z[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
