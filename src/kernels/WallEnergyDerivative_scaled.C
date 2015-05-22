/**
 * @file   WallEnergyDerivative.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 11:59:56 2013
 *
 * @brief
 *
 *
 */

#include "WallEnergyDerivative_scaled.h"

template<>
InputParameters validParams<WallEnergyDerivative_scaled>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1"," "); //FIXME: Give me an explanation
  params.addRequiredParam<Real>("G110"," "); //FIXME: Give me an explanation
  params.addRequiredParam<Real>("G11/G110"," ");
  params.addRequiredParam<Real>("G12/G110"," ");
  params.addRequiredParam<Real>("G44/G110"," ");
  params.addRequiredParam<Real>("G44P/G110"," ");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  params.addParam<Real>("energy_scale",1.0,"energy scale");
  return params;
}


//Constructor
WallEnergyDerivative_scaled::WallEnergyDerivative_scaled(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_i_grad((_component==0)? coupledGradient("polar_x") :(_component==1)? coupledGradient("polar_y"): coupledGradient("polar_z")),
   _polar_j_grad((_component==0)? coupledGradient("polar_y"): (_component==1)? coupledGradient("polar_z"): coupledGradient("polar_x")),
   _polar_k_grad((_component==0)? coupledGradient("polar_z"): (_component==1)? coupledGradient("polar_x"): coupledGradient("polar_y")),
   _ii(_component),
   _jj((_component==0)? 1 : (_component==1)? 2: 0),
   _kk((_component==0)? 2 : (_component==1)? 0: 1),
   _alpha1(getParam<Real>("alpha1")),
   _G110(getParam<Real>("G110")),
   _G11(getParam<Real>("G11/G110")*_G110),
   _G12(getParam<Real>("G12/G110")*_G110),
   _G44(getParam<Real>("G44/G110")*_G110),
   _G44P(getParam<Real>("G44P/G110")*_G110),
   _len_scale(getParam<Real>("len_scale")),
   _energy_scale(getParam<Real>("energy_scale"))
{
  //only for debug purpose
  std::cout<<"_G110="<<_G110<<"\n";
  std::cout<<"_G11="<<_G11<<"\n";
  std::cout<<"_G12="<<_G12<<"\n";
  std::cout<<"_G44="<<_G44<<"\n";
  std::cout<<"_G44P="<<_G44P<<"\n";
}


//TODO:Overload functions

Real
WallEnergyDerivative_scaled::computeQpResidual()
{
  return pow(_len_scale,2)*(_G11*_polar_i_grad[_qp](_ii)*_grad_test[_i][_qp](_ii)+
    _G12*(_polar_j_grad[_qp](_jj)+_polar_k_grad[_qp](_kk))*_grad_test[_i][_qp](_ii)+
    _G44*(_polar_i_grad[_qp](_jj)+_polar_j_grad[_qp](_ii))*_grad_test[_i][_qp](_jj)+ 
    _G44*(_polar_i_grad[_qp](_kk)+_polar_k_grad[_qp](_ii))*_grad_test[_i][_qp](_kk)+
    _G44P*(_polar_i_grad[_qp](_jj)-_polar_j_grad[_qp](_ii))*_grad_test[_i][_qp](_jj)+
    _G44P*(_polar_i_grad[_qp](_kk)-_polar_k_grad[_qp](_ii))*_grad_test[_i][_qp](_kk));

}

Real
WallEnergyDerivative_scaled::computeQpJacobian()
{
  return pow(_len_scale,2)*(_G11*_grad_phi[_j][_qp](_ii)*_grad_test[_i][_qp](_ii)+
          (_G44+_G44P)*_grad_phi[_j][_qp](_jj)*_grad_test[_i][_qp](_jj)+
	  (_G44+_G44P)*_grad_phi[_j][_qp](_kk)*_grad_test[_i][_qp](_kk));
}

Real
WallEnergyDerivative_scaled::computeQpOffDiagJacobian(unsigned int jvar)
{
  mooseAssert(jvar!=variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var)
  {
    const unsigned int _jj = (jvar==_polar_x_var)? 0: (jvar==_polar_y_var)? 1 : 2;
    return pow(_len_scale,2)*(_G12*_grad_phi[_j][_qp](_jj)*_grad_test[_i][_qp](_ii)+(_G44-_G44P)*_grad_phi[_j][_qp](_ii)*_grad_test[_i][_qp](_jj));
  }else{
    return 0.0;
  }
}
