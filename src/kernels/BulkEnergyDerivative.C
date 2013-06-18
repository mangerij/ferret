/**
 * @file   BulkEnergyDerivative.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 11:59:56 2013
 * 
 * @brief  
 * 
 * 
 */

#include "BulkEnergyDerivative.h"
#include<cmath>

template<>
InputParameters validParams<BulkEnergyDerivative>()
{
  InputParameters params = validParams<Kernel>();
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
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}



//Constructor
BulkEnergyDerivative::BulkEnergyDerivative(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _alpha1(getParam<Real>("alpha1")),
   _alpha11(getParam<Real>("alpha11")),
   _alpha12(getParam<Real>("alpha12")),
   _alpha111(getParam<Real>("alpha111")),
   _alpha112(getParam<Real>("alpha112")),
   _alpha123(getParam<Real>("alpha123"))
   
{
  
} 


//TODO:Overload functions

Real
BulkEnergyDerivative::computeQpResidual()
{
  const VariableValue& _polar_i= (_component==0)? _polar_x : (_component==1)? _polar_y: _polar_z;
  const VariableValue& _polar_j= (_component==0)? _polar_y : (_component==1)? _polar_z: _polar_x;
  const VariableValue& _polar_k= (_component==0)? _polar_z : (_component==1)? _polar_x: _polar_y;
  return (2*_alpha1*_polar_i[_qp]+
	  4*_alpha11*pow(_polar_i[_qp],3)+
	  2*_alpha12*_polar_i[_qp]*(pow(_polar_j[_qp],2)+pow(_polar_k[_qp],2))+
	  6*_alpha111*pow(_polar_i[_qp],5)+
	  4*_alpha112*pow(_polar_i[_qp],3)*(_polar_j[_qp]*_polar_j[_qp]+_polar_k[_qp]*_polar_k[_qp])+
	  2*_alpha112*_polar_i[_qp]*(pow(_polar_j[_qp],4)+pow(_polar_k[_qp],3))+
	  2*_alpha123*_polar_i[_qp]*pow(_polar_j[_qp],2)*pow(_polar_k[_qp],2))*_test[_i][_qp];
  
}

Real
BulkEnergyDerivative::computeQpJacobian()
{
  const VariableValue& _polar_i= (_component==0)? _polar_x : (_component==1)? _polar_y: _polar_z;
  const VariableValue& _polar_j= (_component==0)? _polar_y : (_component==1)? _polar_z: _polar_x;
  const VariableValue& _polar_k= (_component==0)? _polar_z : (_component==1)? _polar_x: _polar_y;
  return (2*_alpha1+12*_alpha11*pow(_polar_i[_qp],2)+
	  2*_alpha12*(pow(_polar_j[_qp],2)+pow(_polar_k[_qp],2))+30*_alpha111*pow(_polar_i[_qp],4)+
	  12*_alpha112*pow(_polar_i[_qp],2)*(pow(_polar_j[_qp],2)+pow(_polar_k[_qp],2))+2*_alpha112*(pow(_polar_j[_qp],4)+pow(_polar_k[_qp],4))+
	  2*_alpha123*pow(_polar_j[_qp],2)*pow(_polar_k[_qp],2)
	  )*_test[_i][_qp]*_phi[_j][_qp];
}

Real
BulkEnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  mooseAssert(jvar!=variable().index(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var){
    const VariableValue& _polar_i= (_component==0)? _polar_x : (_component==1)? _polar_y: _polar_z;
    const VariableValue& _polar_j= (jvar==_polar_x_var)? _polar_x : (jvar==_polar_y_var)? _polar_y: _polar_z;
    const VariableValue& _polar_k= ((_component==0 && jvar==_polar_y_var) || (_component==1 && jvar==_polar_x_var) )? _polar_z : ( (_component==0 && jvar==_polar_z_var) || (_component==2 && jvar==_polar_x_var))? _polar_y: _polar_z;
    return (4*_alpha12*_polar_i[_qp]*_polar_j[_qp]
	  +8*_alpha112*pow(_polar_i[_qp],3)*_polar_j[_qp]+8*_alpha112*_polar_i[_qp]*pow(_polar_j[_qp],3)
	  +4*_alpha123*_polar_i[_qp]*_polar_j[_qp]*pow(_polar_k[_qp],2)
	  )*_test[_i][_qp]*_phi[_j][_qp];
  }else
    return 0.0;
}
