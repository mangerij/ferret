/**
 * @file   TotalEnergyGradient.h
 * @author S. Gu <mangerij@anl.gov>
 * @brief
 *
 *
 */


#ifndef TOTALENERGYGRADIENT_H
#define TOTALENERGYGRADIENT_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class TotalEnergyGradient;

template<>
InputParameters validParams<TotalEnergyGradient>();

//TODO: change the base class!
class TotalEnergyGradient : public ElementIntegralPostprocessor
{
public:
  TotalEnergyGradient(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpIntegral();
  const unsigned int _component;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableGradient& _polar_i_grad;
  const VariableGradient& _polar_j_grad;
  const VariableGradient& _polar_k_grad;
  const unsigned int _ii, _jj, _kk;
  const Real _G110,_G11, _G12, _G44, _G44P;
  const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;
  const Real _len_scale;
  const Real _energy_scale;
};

#endif
