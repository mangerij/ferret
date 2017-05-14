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

#include "ComputeElasticityTensor.h"
#include "ComputeRotatedElasticityTensorBase.h"
#include "ComputePhotostrictiveTensor.h"
#include "RotationTensor.h"
#include "RankFourTensor.h"

template<>
InputParameters validParams<ComputePhotostrictiveTensor>()
{
  InputParameters params = validParams<ComputeRotatedPhotostrictiveTensorBase>();
  params.addClassDescription("Compute a photostrictive tensor.");
  params.addRequiredParam<std::vector<Real> >("P_mnkl", "elasto-optic tensor for material");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputePhotostrictiveTensor::ComputePhotostrictiveTensor(const InputParameters & parameters) :
    ComputeRotatedPhotostrictiveTensorBase(parameters),
    _Pmnkl(getParam<std::vector<Real> >("P_mnkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _Pmnkl.rotate(R);
}

void
ComputePhotostrictiveTensor::computeQpPhotostrictiveTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _photostrictive_tensor[_qp] = _Pmnkl;
}


//void
//ComputePhotostrictiveTensor::computeQpUnstrainedRefractiveIndex()
//{
//  // Assume that n_e is along the z-axis for now
//  // note the regular birefringence is quantified by _ne - _no
//  RealVectorValue n(_no, _no, _ne); 

//  // Rotate the indicatrix such that it is aligned with the crystallographic direction A_i = R_{ij} A_j
//  RealVectorValue nR(R(0, 0) * n(0) + R(0, 1) * n(1) + R(0, 2) * n(2), R(1, 0) * n(0) + R(1, 1) * n(1) + R(1, 2) * n(2), R(2, 0) * n(0) + R(2, 1) * n(1) + R(2, 2) * n(2));1
//}

//void
//ComputePhotostrictiveTensor::computeQpStrainedRefractiveIndex()
//{
//  // Assume that n_e is along the z-axis for now
//  // note the regular birefringence is quantified by _ne - _no
//  RealVectorValue n(_no, _no, _ne); 
//
//  // Rotate the indicatrix such that it is aligned with the crystallographic direction A_i = R_{ij} A_j
//  RealVectorValue nR(R(0, 0) * n(0) + R(0, 1) * n(1) + R(0, 2) * n(2), R(1, 0) * n(0) + R(1, 1) * n(1) + R(1, 2) * n(2), R(2, 0) * n(0) + R(2, 1) * n(1) + R(2, 2) * n(2));1
//
//
//  then store in an aux kernel 
//}
