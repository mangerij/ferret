/**
 * @file   SinFuncAux.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Aug 12 12:15:52 2013
 *
 * @brief
 *
 *
 */


#include "SinFuncAux.h"

template<>
InputParameters validParams<SinFuncAux>()
{
  InputParameters params = validParams<AuxKernel>();
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

SinFuncAux::SinFuncAux(const std::string & name, InputParameters parameters) :
  AuxKernel(name, parameters),
  _amplitude(getParam<Real>("amplitude")),
  _wave_length_x(getParam<Real>("wave_length_x")),
  _wave_length_y(getParam<Real>("wave_length_y")),
  _wave_length_z(getParam<Real>("wave_length_z")),
  _phrase_x(getParam<Real>("phrase_x")),
  _phrase_y(getParam<Real>("phrase_y")),
  _phrase_z(getParam<Real>("phrase_z")),
  _vertical_shift(getParam<Real>("vertical_shift"))
{}

/** TODO: Overload functions:
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real
SinFuncAux::computeValue()
{
   Point p=_q_point[_qp];
   Real rv=_vertical_shift+_amplitude*sin((2*pi)/_wave_length_x*p(0)+_phrase_x)*sin((2*pi)/_wave_length_y*p(1)+_phrase_y)*sin((2*pi)/_wave_length_z*p(2)+_phrase_z);
  return rv;
}
