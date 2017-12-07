#include "DomainFunc.h"
#include "MooseRandom.h"

template <>
InputParameters
validParams<DomainFunc>()
{
  InputParameters params = validParams<Function>();
  params.addRequiredParam<Real>("ax", "Width of ferroelectric");
  params.addRequiredParam<Real>("af", "Thickness of ferroelectric");
  params.addRequiredParam<Real>("min", "Min value for random perturbation of polarization");
  params.addRequiredParam<Real>("max", "Max value for random perturbation of polarization");
  return params;
}

DomainFunc::DomainFunc(const InputParameters & parameters)
  : Function(parameters), 
  _ax(getParam<Real>("ax")),
  _af(getParam<Real>("af")),
  _min(getParam<Real>("min")),
  _max(getParam<Real>("max"))
{
    MooseRandom::seed(0.0);
}

Real
DomainFunc::value(Real /*t*/, const Point & p)
{
  Real rand_num = MooseRandom::rand();

  // Between 0 and range
  rand_num *= (_max - _min);

  // Between min and max
  rand_num += _min;
    
  return 1.4 * std::cos(2.0 * libMesh::pi * p(0) / _ax) * std::cos(libMesh::pi * p(2) / _af) + rand_num;
}
