/**
 * @file   PerturbedIC.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Jun 26 18:39:40 2013
 *
 * @brief Perturbed a given number randomly.
 *
 *
 */
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



PerturbedIC::PerturbedIC(const std::string & name, InputParameters parameters) :
  InitialCondition(name, parameters),
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
