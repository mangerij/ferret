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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "PerturbedIC.h"
#include "MooseRandom.h"
#include "libmesh/point.h"

template<>
InputParameters validParams<PerturbedIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("mean", "the mean value");
  params.addRequiredParam<Real>("factor", "multiplicative perturbation factor");
  params.addParam<unsigned int>("seed",0,"seed for the random number generator");
  return params;
}



PerturbedIC::PerturbedIC(const InputParameters & parameters) :
  InitialCondition(parameters),
  _mean(getParam<Real>("mean")),
  _factor(getParam<Real>("factor"))
{
  MooseRandom::seed(getParam<unsigned int>("seed"));
}

Real
PerturbedIC::value(const Point & /*p*/)
{
  Real rand_num=MooseRandom::rand();
  return _mean+2*(rand_num-0.5)*_factor;
}
