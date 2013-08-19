/**
 * @file   RandomFunc.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Aug 12 16:29:47 2013
 *
 * @brief
 *
 *
 */

#include "RandomFunc.h"
#include "MooseRandom.h"

template<>
InputParameters validParams<RandomFunc>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("min", 0.0, "Lower bound of the randomly generated values");
  params.addParam<Real>("max", 1.0, "Upper bound of the randomly generated values");
  params.addParam<unsigned int>("seed", 0, "Seed value for the random number generator");
  return params;
}

RandomFunc::RandomFunc(const std::string & name, InputParameters parameters) :
  Function(name, parameters),
  _min(getParam<Real>("min")),
  _max(getParam<Real>("max")),
  _range(_max - _min)
{}

Real
  RandomFunc::value(Real /*t*/, const Point & /*p*/)
{
   //Random number between 0 and 1
  Real rand_num = MooseRandom::rand();

  //Between 0 and range
  rand_num *= _range;

  //Between min and max
  rand_num += _min;

  return rand_num;
}
