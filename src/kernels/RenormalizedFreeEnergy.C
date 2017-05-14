/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "RenormalizedFreeEnergy.h"
#include<cmath>

template<>
InputParameters validParams<RenormalizedFreeEnergy>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<Real>("epsilon", "Constant value for strain");
  params.addRequiredParam<Real>("T", "the temperature");
  return params;
}

RenormalizedFreeEnergy::RenormalizedFreeEnergy(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _epsilon(getParam<Real>("epsilon")),
   _T(getParam<Real>("T"))
{
  std::cout<<"_epsilon ="<<_epsilon<<"\n";
}

Real
RenormalizedFreeEnergy::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (3.08064 * std::pow(_polar_x[_qp], 3.0) + 1.56 * std::pow(_polar_x[_qp], 5.0) + 0.977049 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0) + 0.507273 * _polar_x[_qp] * _polar_z[_qp] * _polar_z[_qp] - 7.4 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 0.61 * (4 * std::pow(_polar_x[_qp], 3.0) * std::pow(_polar_y[_qp], 2.0) + 4 * std::pow(_polar_x[_qp], 3.0) * std::pow(_polar_z[_qp], 2.0) + 2 * _polar_x[_qp] * (std::pow(_polar_y[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0))) + 2 * _polar_x[_qp] * (-4.92223 * (-765.1 + _T) * std::pow(10.0, -7.0) - 19.0909 * _epsilon));
   //3.08064 Px^3 + 1.56 Px^5 + 0.977049 Px Py^2 + 0.507273 Px Pz^2 - 7.4 Px Py^2 Pz^2 + 0.61 (4 Px^3 Py^2 + 4 Px^3 Pz^2 + 2 Px (Py^4 + Pz^4)) + 2 Px (0.000229917 - 19.0909 \[Epsilon])
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (3.08064 * std::pow(_polar_y[_qp], 3.0) + 1.56 * std::pow(_polar_y[_qp], 5.0) + 0.977049 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0) + 0.507273 * _polar_y[_qp] * _polar_z[_qp] * _polar_z[_qp] - 7.4 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 0.61 * (4 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 3.0) + 4 * std::pow(_polar_y[_qp], 3.0) * std::pow(_polar_z[_qp], 2.0) + 2 * _polar_y[_qp] * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0))) + 2 * _polar_y[_qp] * (-4.92223 * (-765.1 + _T) * std::pow(10.0, -7.0) - 19.0909 * _epsilon));
    //0.977049 Px^2 Py + 3.08064 Py^3 + 1.56 Py^5 + 0.507273 Py Pz^2 - 7.4 Px^2 Py Pz^2 + 0.61 (4 Px^2 Py^3 + 4 Py^3 Pz^2 + 2 Py (Px^4 + Pz^4)) +  2 Py (0.000229917 - 19.0909 \[Epsilon])
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (-7.4 * _polar_z[_qp] * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + 0.529274 * std::pow(_polar_z[_qp], 3.0) + 1.56 * std::pow(_polar_z[_qp], 5.0) + 0.253636 * (2 * std::pow(_polar_x[_qp], 2.0) * _polar_z[_qp] + 2 * std::pow(_polar_y[_qp], 2.0) * _polar_z[_qp]) + 0.61 * (2 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) * _polar_z[_qp] + 4 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 3.0) + 4 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 3.0)) + 2 * _polar_z[_qp] * (-4.92223 * (-765.1 + _T) * std::pow(10.0, -7.0) + 15.7576 * _epsilon));
    //-7.4 Px^2 Py^2 Pz + 0.529274 Pz^3 + 1.56 Pz^5 +  0.253636 (2 Px^2 Pz + 2 Py^2 Pz) +  0.61 (2 (Px^4 + Py^4) Pz + 4 Px^2 Pz^3 + 4 Py^2 Pz^3) +  2 Pz (0.000229917 + 15.7576 \[Epsilon])
  }
  else 
    return 0.0;
}

Real
RenormalizedFreeEnergy::computeQpJacobian()
{
  if (_component == 0)
  {
    return _phi[_j][_qp] * _test[_i][_qp] * (9.24192 * std::pow(_polar_x[_qp], 2.0) + 7.8 * std::pow(_polar_x[_qp], 4.0) + 0.977049 * std::pow(_polar_y[_qp], 2.0) + 0.507273 * std::pow(_polar_z[_qp], 2.0) - 7.4 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 0.61 * ( 12 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + 12.0 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 2.0 * std::pow(_polar_y[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0)) + 2 * (-4.92223 * (-765.1 + _T) * std::pow(10.0, -7.0) - 19.0909 * _epsilon));
  //9.24192 Px^2 + 7.8 Px^4 + 0.977049 Py^2 + 0.507273 Pz^2 -  7.4 Py^2 Pz^2 +  0.61 (12 Px^2 Py^2 + 12 Px^2 Pz^2 + 2 (Py^4 + Pz^4)) +  2 (0.000229917 - 19.0909 \[Epsilon])
  }
  else if (_component == 1)
  {
    return _phi[_j][_qp] * _test[_i][_qp] * (9.24192 * std::pow(_polar_y[_qp], 2.0) + 7.8 * std::pow(_polar_y[_qp], 4.0) + 0.977049 * std::pow(_polar_x[_qp], 2.0) + 0.507273 * std::pow(_polar_z[_qp], 2.0) - 7.4 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 0.61 * ( 12 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_x[_qp], 2.0) + 12.0 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 2.0 * std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0)) + 2 * (-4.92223 * (-765.1 + _T) * std::pow(10.0, -7.0) - 19.0909 * _epsilon));
  //0.977049 Px^2 + 9.24192 Py^2 + 7.8 Py^4 + 0.507273 Pz^2 -  7.4 Px^2 Pz^2 +  0.61 (12 Px^2 Py^2 + 12 Py^2 Pz^2 + 2 (Px^4 + Pz^4)) +  2 (0.000229917 - 19.0909 \[Epsilon])
  }
  else if (_component == 2)
  {
    return _phi[_j][_qp] * _test[_i][_qp] * (-7.4 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_x[_qp], 2.0) + 0.253636 * ( 2 * std::pow(_polar_x[_qp], 2.0) + 2 * std::pow(_polar_y[_qp], 2.0) ) + 1.58782 * std::pow(_polar_z[_qp], 2.0) + 7.8 * std::pow(_polar_z[_qp], 4.0) + 0.61 *  (2 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) + 12 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 12 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0)) + 2 * (-4.92223 * (-765.1 + _T) * std::pow(10.0, -7.0) + 15.7576 * _epsilon));
  //-7.4 Px^2 Py^2 + 0.253636 (2 Px^2 + 2 Py^2) + 1.58782 Pz^2 +  7.8 Pz^4 + 0.61 (2 (Px^4 + Py^4) + 12 Px^2 Pz^2 + 12 Py^2 Pz^2) +  2 (0.000229917 + 15.7576 \[Epsilon])
  }
  else 
    return 0.0;
}

Real
RenormalizedFreeEnergy::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * _test[_i][_qp] * (1.9541 * _polar_x[_qp] * _polar_y[_qp] + 0.61 * (8.0 * std::pow(_polar_x[_qp], 3.0) * _polar_y[_qp] + 8 * std::pow(_polar_y[_qp], 3.0) * _polar_x[_qp]) - 14.8 * _polar_x[_qp] * _polar_y[_qp] * std::pow(_polar_z[_qp], 2.0));
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * _test[_i][_qp] * (1.01455 * _polar_x[_qp] * _polar_z[_qp] + 0.61 * (8.0 * std::pow(_polar_x[_qp], 3.0) * _polar_z[_qp] + 8 * std::pow(_polar_z[_qp], 3.0) * _polar_x[_qp]) - 14.8 * _polar_x[_qp] * _polar_z[_qp] * std::pow(_polar_y[_qp], 2.0));
    }
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * _test[_i][_qp] * (1.9541 * _polar_x[_qp] * _polar_y[_qp] + 0.61 * (8.0 * std::pow(_polar_x[_qp], 3.0) * _polar_y[_qp] + 8 * std::pow(_polar_y[_qp], 3.0) * _polar_x[_qp]) - 14.8 * _polar_x[_qp] * _polar_y[_qp] * std::pow(_polar_z[_qp], 2.0));
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * _test[_i][_qp] * (1.01455 * _polar_y[_qp] * _polar_z[_qp] + 0.61 * (8.0 * std::pow(_polar_y[_qp], 3.0) * _polar_z[_qp] + 8 * std::pow(_polar_z[_qp], 3.0) * _polar_y[_qp]) - 14.8 * _polar_y[_qp] * _polar_z[_qp] * std::pow(_polar_x[_qp], 2.0));
    }
    else
      return 0.0;
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * _test[_i][_qp] * (1.01455 * _polar_x[_qp] * _polar_z[_qp] + 0.61 * (8.0 * std::pow(_polar_x[_qp], 3.0) * _polar_z[_qp] + 8 * std::pow(_polar_z[_qp], 3.0) * _polar_x[_qp]) - 14.8 * _polar_x[_qp] * _polar_z[_qp] * std::pow(_polar_y[_qp], 2.0));
    }
    else if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * _test[_i][_qp] * (1.01455 * _polar_y[_qp] * _polar_z[_qp] + 0.61 * (8.0 * std::pow(_polar_y[_qp], 3.0) * _polar_z[_qp] + 8 * std::pow(_polar_z[_qp], 3.0) * _polar_y[_qp]) - 14.8 * _polar_y[_qp] * _polar_z[_qp] * std::pow(_polar_x[_qp], 2.0));
    }
    else
      return 0.0;
  }
  else 
    return 0.0;
}




