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

#include "BulkEnergyDerivativeSixth.h"
#include "libmesh/utility.h"
#include <vector>

registerMooseObject("FerretApp", BulkEnergyDerivativeSixth);

template<>
InputParameters validParams<BulkEnergyDerivativeSixth>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for the local free energy which is an sixth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<std::vector<Real>>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<std::vector<Real>>("alpha11", "The coefficients of the Landau expansion");
  params.addRequiredParam<std::vector<Real>>("alpha12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

BulkEnergyDerivativeSixth::BulkEnergyDerivativeSixth(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _alpha1(getParam<std::vector<Real>>("alpha1").at(0)),
   _alpha3(getParam<std::vector<Real>>("alpha1").size() > 1 ? getParam<std::vector<Real>>("alpha1").at(1) : _alpha1),
   _alpha11(getParam<std::vector<Real>>("alpha11").at(0)),
   _alpha33(getParam<std::vector<Real>>("alpha11").size() > 1 ? getParam<std::vector<Real>>("alpha11").at(1) : _alpha11),
   _alpha12(getParam<std::vector<Real>>("alpha12").at(0)),
   _alpha13(getParam<std::vector<Real>>("alpha12").size() > 1 ? getParam<std::vector<Real>>("alpha12").at(1) : _alpha12),
   _alpha111(getParam<Real>("alpha111")),
   _alpha112(getParam<Real>("alpha112")),
   _alpha123(getParam<Real>("alpha123")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkEnergyDerivativeSixth::computeQpResidual()
{
  if (_component == 0)
  {
    return (
        2.0 * _alpha1 * _polar_x[_qp] +
        4.0 * _alpha11 * Utility::pow<3>(_polar_x[_qp]) +
        2.0 * _alpha13 * _polar_x[_qp] * Utility::pow<2>(_polar_z[_qp]) +
        2.0 * _alpha12 * _polar_x[_qp] * Utility::pow<2>(_polar_y[_qp]) +
        6.0 * _alpha111 * Utility::pow<5>(_polar_x[_qp]) +
        _alpha112 * (
            4.0 * Utility::pow<3>(_polar_x[_qp]) * ( Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) ) +
            2.0 * _polar_x[_qp] * ( Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp]) )
        ) +
        2.0 * _alpha123 * _polar_x[_qp] * Utility::pow<2>(_polar_y[_qp]) * Utility::pow<2>(_polar_z[_qp])
    )  * _test[_i][_qp];
  }
  else if (_component == 1)
  {
    return (
        2.0 * _alpha1 * _polar_y[_qp] +
        4.0 * _alpha11 * Utility::pow<3>(_polar_y[_qp]) +
        2.0 * _alpha13 * _polar_y[_qp] * Utility::pow<2>(_polar_z[_qp]) +
        2.0 * _alpha12 * _polar_y[_qp] * Utility::pow<2>(_polar_x[_qp]) +
        6.0 * _alpha111 * Utility::pow<5>(_polar_y[_qp]) +
        _alpha112 * (
            4.0 * Utility::pow<3>(_polar_y[_qp]) * ( Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]) ) +
            2.0 * _polar_y[_qp] * ( Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_z[_qp]) )
        ) +
        2.0 * _alpha123 * _polar_y[_qp] * Utility::pow<2>(_polar_x[_qp]) * Utility::pow<2>(_polar_z[_qp])
    )  * _test[_i][_qp];
  }
  else if (_component == 2)
  {
    return (
        2.0 * _alpha3 * _polar_z[_qp] +
        4.0 * _alpha33 * Utility::pow<3>(_polar_z[_qp]) +
        2.0 * _alpha13 * _polar_z[_qp] * ( Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) ) +
        6.0 * _alpha111 * Utility::pow<5>(_polar_z[_qp]) +
        _alpha112 * (
            4.0 * Utility::pow<3>(_polar_z[_qp]) * ( Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) ) +
            2.0 * _polar_z[_qp] * ( Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) ) 
        ) +
        2.0 * _alpha123 * _polar_z[_qp] * Utility::pow<2>(_polar_x[_qp]) * Utility::pow<2>(_polar_y[_qp])
    )  * _test[_i][_qp];
  }
  else 
  {
    return 0.0;
  }
}

Real
BulkEnergyDerivativeSixth::computeQpJacobian()
{  
  if (_component == 0)
  {
    return (
        2.0 * _alpha1 + 12.0 * _alpha11 * Utility::pow<2>(_polar_x[_qp]) + 2.0 * _alpha13 * Utility::pow<2>(_polar_z[_qp]) +
        2.0 * _alpha12 * Utility::pow<2>(_polar_y[_qp]) + 30.0 * _alpha111 * Utility::pow<4>(_polar_x[_qp]) +
        _alpha112 * ( 
            12.0 * Utility::pow<2>(_polar_x[_qp]) * ( Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) ) + 
            2.0 * ( Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp]) ) 
        ) + 
        2.0 * _alpha123 * Utility::pow<2>(_polar_y[_qp]) * Utility::pow<2>(_polar_z[_qp])
    )  * _phi[_j][_qp] * _test[_i][_qp];
  }
  else if (_component == 1)
  {
    return (
        2.0 * _alpha1 + 12.0 * _alpha11 * Utility::pow<2>(_polar_y[_qp]) + 2.0 * _alpha13 * Utility::pow<2>(_polar_z[_qp]) +
        2.0 * _alpha12 * Utility::pow<2>(_polar_x[_qp]) + 30.0 * _alpha111 * Utility::pow<4>(_polar_y[_qp]) +
        _alpha112 * ( 
            12.0 * Utility::pow<2>(_polar_y[_qp]) * ( Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]) ) + 
            2.0 * ( Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_z[_qp]) ) 
        ) + 
        2.0 * _alpha123 * Utility::pow<2>(_polar_x[_qp]) * Utility::pow<2>(_polar_z[_qp])
    )  * _phi[_j][_qp] * _test[_i][_qp];
  }
  else if (_component == 2)
  {
    return (
        2.0 * _alpha3 + 12.0 * _alpha33 * Utility::pow<2>(_polar_z[_qp]) + 
        2.0 * _alpha13 * ( Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) ) + 
        30.0 * _alpha111 * Utility::pow<4>(_polar_z[_qp]) +
        _alpha112 * ( 
            12.0 * Utility::pow<2>(_polar_z[_qp]) * ( Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) ) + 
            2.0 * ( Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) ) 
        ) + 
        2.0 * _alpha123 * Utility::pow<2>(_polar_x[_qp]) * Utility::pow<2>(_polar_y[_qp])
    )  * _phi[_j][_qp] * _test[_i][_qp];
  }
  else 
  {
    return 0.0;
  }
}


Real
BulkEnergyDerivativeSixth::computeQpOffDiagJacobian(unsigned int jvar)
{
  if(_component == 0) {
    if (jvar == _polar_y_var)
    {
        return (
            4.0 * _alpha12 * _polar_x[_qp] * _polar_y[_qp] + 
            _alpha112 * ( 
                8.0 * Utility::pow<3>(_polar_x[_qp]) * _polar_y[_qp] + 
                8.0 * _polar_x[_qp] *  Utility::pow<3>(_polar_y[_qp]) 
            ) + 
            4.0 * _alpha123 * _polar_x[_qp] * _polar_y[_qp] * Utility::pow<2>(_polar_z[_qp])
        ) * _phi[_j][_qp] * _test[_i][_qp];
    }
    else if (jvar == _polar_z_var)
    {
        return (
            4.0 * _alpha13 * _polar_x[_qp] * _polar_z[_qp] +
            _alpha112 * ( 
                8.0 * Utility::pow<3>(_polar_x[_qp]) * _polar_z[_qp] + 
                8.0 * _polar_x[_qp] * Utility::pow<3>(_polar_z[_qp]) 
            ) + 
            4.0 * _alpha123 * _polar_x[_qp] * Utility::pow<2>(_polar_y[_qp]) * _polar_z[_qp]
        ) * _phi[_j][_qp] * _test[_i][_qp];
    }
    else
    {
      return 0.0;
    }
  }
  else if(_component == 1) {
    if(jvar == _polar_x_var) {
        return (
            4.0 * _alpha12 * _polar_y[_qp] * _polar_x[_qp] +
            _alpha112 * ( 
                8.0 * Utility::pow<3>(_polar_y[_qp]) * _polar_x[_qp] + 
                8.0 * _polar_y[_qp] * Utility::pow<3>(_polar_x[_qp])
            ) + 
            4.0 * _alpha123 * _polar_y[_qp] * _polar_x[_qp] * Utility::pow<2>(_polar_z[_qp])
        ) * _phi[_j][_qp] * _test[_i][_qp];
    } 
    else if(jvar == _polar_z_var) {
        return (
            4.0 * _alpha13 * _polar_y[_qp] * _polar_z[_qp] +
            _alpha112 * ( 
                8.0 * Utility::pow<3>(_polar_y[_qp]) * _polar_z[_qp] + 
                8.0 * _polar_y[_qp] * Utility::pow<3>(_polar_z[_qp]) 
            ) + 
            4.0 * _alpha123 * _polar_y[_qp] * Utility::pow<2>(_polar_x[_qp]) * _polar_z[_qp]
        ) * _phi[_j][_qp] * _test[_i][_qp];
    }
    else 
    {
        return 0.0;
    }
  }
  else if (_component == 2) {
    if(jvar == _polar_x_var) {
        return (
            4.0 * _alpha13 * _polar_z[_qp] * _polar_x[_qp] + 
            _alpha112 * ( 
                8.0 * Utility::pow<3>(_polar_z[_qp]) * _polar_x[_qp] + 
                8.0 * _polar_z[_qp] * Utility::pow<3>(_polar_x[_qp])
            ) + 
            4.0 * _alpha123 * _polar_z[_qp] * _polar_x[_qp] * Utility::pow<2>(_polar_y[_qp])
        ) * _phi[_j][_qp] * _test[_i][_qp];
    } 
    else if(jvar == _polar_y_var) {
        return (
            4.0 * _alpha13 * _polar_z[_qp] * _polar_y[_qp] + 
            _alpha112 * ( 
                8.0 * Utility::pow<3>(_polar_z[_qp]) * _polar_y[_qp] + 
                8.0 * _polar_z[_qp] * Utility::pow<3>(_polar_y[_qp])  
            ) + 
            4.0 * _alpha123 * _polar_z[_qp] * Utility::pow<2>(_polar_x[_qp]) * _polar_y[_qp]
        ) * _phi[_j][_qp] * _test[_i][_qp];
    }
    else
    {
        return 0.0;
    }
  }
  else 
  {
    return 0.0;
  }
}
