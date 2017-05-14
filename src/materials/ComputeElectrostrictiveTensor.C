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
#include "ComputeElectrostrictiveTensor.h"
#include "RotationTensor.h"
#include "RankFourTensor.h"

template<>
InputParameters validParams<ComputeElectrostrictiveTensor>()
{
  InputParameters params = validParams<ComputeRotatedElectrostrictiveTensorBase>();
  params.addClassDescription("Compute an electrostrictive tensor.");
  params.addParam<bool>("compute_electrostrictive_coeff", false, "compute the electrostrictive coefficients Q_mnkl");
  params.addRequiredParam<std::vector<Real> >("Q_mnkl", "electrostrictive tensor for material");
  params.addRequiredParam<std::vector<Real> >("C_ijkl", "elastic stiffness tensor for material");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputeElectrostrictiveTensor::ComputeElectrostrictiveTensor(const InputParameters & parameters) :
    ComputeRotatedElectrostrictiveTensorBase(parameters),
    _compute_electrostrictive_coeff(getParam<bool>("compute_electrostrictive_coeff")),
    _Qmnkl(getParam<std::vector<Real> >("Q_mnkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _Cijkl(getParam<std::vector<Real> >("C_ijkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _Qmnkl.rotate(R);
  _Cijkl.rotate(R);
  ///contractions using namespace method
  _qijkl = ElectrostrictiveTensorTools::computeProduct(_Cijkl, _Qmnkl);
}

void
ComputeElectrostrictiveTensor::computeQpElectrostrictiveTensor()
{
  ///Assign an electrostrictive tensor at a given quad point -- in principle we DON'T want this?
  _electrostrictive_tensor[_qp] = _qijkl;
  //Need electrostrictive coefficients stored as well for the polar-optic work. Should be a bool flag here whether we want to store this or not
  if (_compute_electrostrictive_coeff == true)
    _electrostrictive_coefficients[_qp] = _Qmnkl;
}
