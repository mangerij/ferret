/**
 * @file   SinFunc.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Aug 12 16:02:40 2013
 *
 * @brief  compute sinusoidual wave function
 *
 *
 */

#include "SinFunc.h"

template<>
InputParameters validParams<SinFunc>()
{
  InputParameters params = validParams<Function>();
  params.addRequiredParam<Real>("amplitude", "amplitude");
  params.addRequiredParam<Real>("wave_length_x", "wave length in x");
  params.addRequiredParam<Real>("wave_length_y", "wave length in y");
  params.addRequiredParam<Real>("wave_length_z", "wave length in z");
  params.addRequiredParam<Real>("phrase_x","phrase in x");
  params.addRequiredParam<Real>("phrase_y","phrase in y");
  params.addRequiredParam<Real>("phrase_z","phrase in z");
  params.addRequiredParam<Real>("vertical_shift","vertical shift");
  return params;
}

SinFunc::SinFunc(const std::string & name, InputParameters parameters) :
  Function(name, parameters),
  _amplitude(getParam<Real>("amplitude")),
  _wave_length_x(getParam<Real>("wave_length_x")),
  _wave_length_y(getParam<Real>("wave_length_y")),
  _wave_length_z(getParam<Real>("wave_length_z")),
  _phrase_x(getParam<Real>("phrase_x")),
  _phrase_y(getParam<Real>("phrase_y")),
  _phrase_z(getParam<Real>("phrase_z")),
  _vertical_shift(getParam<Real>("vertical_shift"))
{}

Real
SinFunc::value(Real /*t*/, const Point & p)
{
   Real rv=_vertical_shift+_amplitude*sin((2*pi)/_wave_length_x*p(0)+_phrase_x)*sin((2*pi)/_wave_length_y*p(1)+_phrase_y)*sin((2*pi)/_wave_length_z*p(2)+_phrase_z);
   return rv;
}
