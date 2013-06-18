/**
 * @file   PolarMaterial.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun  4 17:17:46 2013
 * 
 * @brief  
 * 
 * 
 */


#include "PolarMaterial.h"

template<>
InputParameters validParams<PolarMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1"," "); //FIXME: Give me an explanation
  params.addRequiredParam<Real>("alpha11"," ");
  params.addRequiredParam<Real>("alpha12"," ");
  params.addRequiredParam<Real>("alpha111"," ");
  params.addRequiredParam<Real>("alpha112"," ");
  params.addRequiredParam<Real>("alpha123"," ");
  params.addRequiredParam<Real>("G110","");
  params.addRequiredParam<Real>("G11/G110"," ");
  params.addRequiredParam<Real>("G12/G110"," ");
  params.addRequiredParam<Real>("G44/G110"," ");
  params.addRequiredParam<Real>("G44P/G110"," ");
  return params;
}

PolarMaterial::PolarMaterial(const std::string & name,
                                 InputParameters parameters)
  :Material(name, parameters),
   _polars(declareProperty<std::vector<Real> >("polars")),
   _polar_grads(declareProperty<std::vector<RealGradient> >("polar_grads")),
   _alpha1(declareProperty<Real>("alpha1")),
   _alpha11(declareProperty<Real>("alpha11")),
   _alpha12(declareProperty<Real>("alpha12")),
   _alpha111(declareProperty<Real>("alpha111")),
   _alpha112(declareProperty<Real>("alpha112")),
   _alpha123(declareProperty<Real>("alpha123")),
   _G11(declareProperty<Real>("G11")),
   _G12(declareProperty<Real>("G12")),
   _G44(declareProperty<Real>("G44")),
  _G44P(declareProperty<Real>("G44P")),
  _alpha1_i(getParam<Real>("alpha1")),
  _alpha11_i(getParam<Real>("alpha11")),
  _alpha12_i(getParam<Real>("alpha12")),
  _alpha111_i(getParam<Real>("alpha111")),
  _alpha112_i(getParam<Real>("alpha112")),
  _alpha123_i(getParam<Real>("alpha123")),
  _G110_i(getParam<Real>("G110")),
  _G11_i(getParam<Real>("G11/G110")*_G110_i),
  _G12_i(getParam<Real>("G12/G110")*_G110_i),
  _G44_i(getParam<Real>("G44/G110")*_G110_i),
  _G44P_i(getParam<Real>("G44P/G110")*_G110_i),
  _polar_x_val(coupledValue("polar_x")),
  _polar_y_val(coupledValue("polar_y")),
  _polar_z_val(coupledValue("polar_z")),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z"))
{
  
}

void
PolarMaterial::computeQpProperties()
{
  std::vector<Real> polar(3);
  polar[0]=_polar_x_val[_qp];
  polar[1]=_polar_y_val[_qp];
  polar[2]=_polar_z_val[_qp];
  _polars[_qp]=polar;
  
  std::vector<RealGradient> polar_grad(3);
  polar_grad[0]=_polar_x_grad[_qp];
  polar_grad[1]=_polar_y_grad[_qp];
  polar_grad[2]=_polar_z_grad[_qp];
  _polar_grads[_qp]=polar_grad;
  
  _alpha1[_qp]=_alpha1_i;
  _alpha11[_qp]=_alpha11_i;
  _alpha12[_qp]=_alpha12_i;
  _alpha111[_qp]=_alpha111_i;
  _alpha112[_qp]=_alpha112_i;
  _alpha123[_qp]=_alpha123_i;
  _G11[_qp]=_G11_i;
  _G12[_qp]=_G12_i;
  _G44[_qp]=_G44_i;
  _G44P[_qp]=_G44P_i;
}
