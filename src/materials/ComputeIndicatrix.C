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

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "ComputeIndicatrix.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

template<>
InputParameters validParams<ComputeIndicatrix>()
{
  InputParameters params = validParams<ComputeRotatedIndicatrixBase>();
  params.addClassDescription("Compute the impermeability tensor.");
  params.addRequiredParam<Real>("n_a", "alpha refractive index");
  params.addRequiredParam<Real>("n_b", "beta refractive index");
  params.addRequiredParam<Real>("n_g", "gamma refractive index");
  return params;
}

ComputeIndicatrix::ComputeIndicatrix(const InputParameters & parameters) :
    ComputeRotatedIndicatrixBase(parameters),
   _na(getParam<Real>("n_a")),
   _nb(getParam<Real>("n_b")),
   _ng(getParam<Real>("n_g"))
{
}

void
ComputeIndicatrix::computeQpIndicatrix()
{
  // Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles);
  //define the principle axis impermeability rotated
  _indicatrix[_qp](0,0) = R(0,0) * R(0,0) / (_na * _na) + R(1,0) * R(1,0) / (_nb * _nb) + R(2,0) * R(2,0) / (_ng * _ng);//
  _indicatrix[_qp](0,1) = R(0,0) * R(0,1) / (_na * _na) + R(1,0) * R(1,1) / (_nb * _nb) + R(2,1) * R(2,0) / (_ng * _ng);//
  _indicatrix[_qp](0,2) = R(0,0) * R(0,2) / (_na * _na) + R(1,0) * R(1,2) / (_nb * _nb) + R(2,0) * R(2,2) / (_ng * _ng);//
  _indicatrix[_qp](1,0) = R(0,1) * R(0,0) / (_na * _na) + R(1,1) * R(1,0) / (_nb * _nb) + R(2,0) * R(2,1) / (_ng * _ng);//
  _indicatrix[_qp](1,1) = R(0,1) * R(0,1) / (_na * _na) + R(1,1) * R(1,1) / (_nb * _nb) + R(2,1) * R(2,1) / (_ng * _ng);//
  _indicatrix[_qp](1,2) = R(0,1) * R(0,2) / (_na * _na) + R(1,1) * R(1,2) / (_nb * _nb) + R(2,1) * R(2,2) / (_ng * _ng);//
  _indicatrix[_qp](2,0) = R(0,2) * R(0,0) / (_na * _na) + R(1,0) * R(1,2) / (_nb * _nb) + R(2,0) * R(2,2) / (_ng * _ng);//
  _indicatrix[_qp](2,1) = R(0,1) * R(0,2) / (_na * _na) + R(1,1) * R(1,2) / (_nb * _nb) + R(2,1) * R(2,2) / (_ng * _ng);//
  _indicatrix[_qp](2,2) = R(0,2) * R(0,2) / (_na * _na) + R(1,2) * R(1,2) / (_nb * _nb) + R(2,2) * R(2,2) / (_ng * _ng);//

  //Moose::out << "\n B"; std::cout << 0; std::cout << 0; Moose::out << " = "; std::cout << _beta_tensor[_qp](0,0);
}


