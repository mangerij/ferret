/**
 * @file   TotalEnergyGradient.C
 * @author J. Mangeri <mangerij@anl.gov>
 *
 * @brief Used to check if the gradient of the energy is indeed at a minimum
 *
 *
 */


#include "TotalEnergyGradient.h"

template<>
InputParameters validParams<TotalEnergyGradient>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1"," "); //FIXME: Give me an explanation
  params.addRequiredParam<Real>("alpha11"," ");
  params.addRequiredParam<Real>("alpha12"," ");
  params.addRequiredParam<Real>("alpha111"," ");
  params.addRequiredParam<Real>("alpha112"," ");
  params.addRequiredParam<Real>("alpha123"," ");
  params.addRequiredParam<Real>("G110"," "); //FIXME: Give me an explanation
  params.addRequiredParam<Real>("G11/G110"," ");
  params.addRequiredParam<Real>("G12/G110"," ");
  params.addRequiredParam<Real>("G44/G110"," ");
  params.addRequiredParam<Real>("G44P/G110"," ");
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  params.addParam<Real>("energy_scale",1.0,"energy scale");
  return params;
}

TotalEnergyGradient::TotalEnergyGradient(const std::string & name, InputParameters parameters) :
  ElementIntegralPostprocessor(name, parameters),
  _component(getParam<unsigned int>("component")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha12(getParam<Real>("alpha12")),
  _alpha111(getParam<Real>("alpha111")),
  _alpha112(getParam<Real>("alpha112")),
  _alpha123(getParam<Real>("alpha123")),
  _polar_x_var(coupled("polar_x")),
  _polar_y_var(coupled("polar_y")),
  _polar_z_var(coupled("polar_z")),
  _polar_i_grad((_component==0)? coupledGradient("polar_x") :(_component==1)? coupledGradient("polar_y"): coupledGradient("polar_z")),
  _polar_j_grad((_component==0)? coupledGradient("polar_y"): (_component==1)? coupledGradient("polar_z"): coupledGradient("polar_x")),
  _polar_k_grad((_component==0)? coupledGradient("polar_z"): (_component==1)? coupledGradient("polar_x"): coupledGradient("polar_y")),
  _ii(_component),
  _jj((_component==0)? 1 : (_component==1)? 2: 0),
  _kk((_component==0)? 2 : (_component==1)? 0: 1),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11/G110")*_G110),
  _G12(getParam<Real>("G12/G110")*_G110),
  _G44(getParam<Real>("G44/G110")*_G110),
  _G44P(getParam<Real>("G44P/G110")*_G110),
  _potential_int_grad(coupledGradient("potential_int")),
  _potential_ext_grad(coupledGradient("potential_ext")),
  _len_scale(getParam<Real>("len_scale")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real

TotalEnergyGradient::computeQpIntegral()
{
  const VariableValue& _polar_i= (_component==0)? _polar_x : (_component==1)? _polar_y: _polar_z;
  const VariableValue& _polar_j= (_component==0)? _polar_y : (_component==1)? _polar_z: _polar_x;
  const VariableValue& _polar_k= (_component==0)? _polar_z : (_component==1)? _polar_x: _polar_y;
  return 
//bulk energy:
((2*_alpha1*_polar_i[_qp]+
	  4*_alpha11*pow(_polar_i[_qp],3)+
	  2*_alpha12*_polar_i[_qp]*(pow(_polar_j[_qp],2)+pow(_polar_k[_qp],2))+
	  6*_alpha111*pow(_polar_i[_qp],5)+
	  4*_alpha112*pow(_polar_i[_qp],3)*(_polar_j[_qp]*_polar_j[_qp]+_polar_k[_qp]*_polar_k[_qp])+
	  2*_alpha112*_polar_i[_qp]*(pow(_polar_j[_qp],4)+pow(_polar_k[_qp],4))+
	   2*_alpha123*_polar_i[_qp]*pow(_polar_j[_qp],2)*pow(_polar_k[_qp],2)))*pow(_len_scale,3.0)*_energy_scale

+(_G11*_polar_i_grad[_qp](_ii)+
    _G12*(_polar_j_grad[_qp](_jj)+_polar_k_grad[_qp](_kk))+
    _G44*(_polar_i_grad[_qp](_jj)+_polar_j_grad[_qp](_ii))+ _G44*(_polar_i_grad[_qp](_kk)+_polar_k_grad[_qp](_ii))+
	   _G44P*(_polar_i_grad[_qp](_jj)-_polar_j_grad[_qp](_ii))+_G44P*(_polar_i_grad[_qp](_kk)-_polar_k_grad[_qp](_ii)))*_len_scale*_energy_scale
//wall energy:

//polar order coupling to electric field
+(0.5*_potential_int_grad[_qp](_component)+_potential_ext_grad[_qp](_component))*pow(_len_scale,2.0)*_energy_scale;
}
