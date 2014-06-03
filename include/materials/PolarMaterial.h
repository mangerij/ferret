/**
 * @file   MaterialTemplate.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun  4 17:02:20 2013
 *
 * @brief
 *
 *
 */

#include "Material.h"

#ifndef POLARMATERIAL_H
#define POLARMATERIAL_H

//Forward Declarations
class PolarMaterial;

template<>
InputParameters validParams<PolarMaterial>();

class PolarMaterial : public Material
{
public:
  PolarMaterial(const std::string & name,
                  InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
  MaterialProperty<std::vector<Real> >& _polars;
  MaterialProperty<std::vector<RealGradient> >& _polar_grads;
  MaterialProperty<Real> & _alpha1, & _alpha11, & _alpha12, & _alpha111, & _alpha112, & _alpha123;
  MaterialProperty<Real> & _G11, & _G12, & _G44, & _G44P;
  const Real _alpha1_i,_alpha11_i, _alpha12_i,_alpha111_i,_alpha112_i,_alpha123_i;
  const Real _G110_i,_G11_i, _G12_i, _G44_i, _G44P_i;
  VariableValue& _polar_x_val;
  VariableValue& _polar_y_val;
  VariableValue& _polar_z_val;
  VariableGradient& _polar_x_grad;
  VariableGradient& _polar_y_grad;
  VariableGradient& _polar_z_grad;
};

#endif //POLARMATERIAL_H
