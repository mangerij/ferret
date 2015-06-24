/**
 * @file   FerroelectricCouplingU.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15. 2015
 * @brief   Implement the kernel for displacement variable corresponding to ferroelectic coupling energy,
 *           Assume the energy has the form -0.5*q_ijkl* ui_j * Pk_l where u is the displacement and P is the polarization.
 */

#include "FerroelectricCouplingU.h"

class FerroelectricCouplingU;

template<>
InputParameters validParams<FerroelectricCouplingU>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}



//Constructor
FerroelectricCouplingU::FerroelectricCouplingU(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _electrostrictive_tensor(getMaterialProperty<ElectrostrictiveTensorR4>("electrostrictive_tensor")),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FerroelectricCouplingU::computeQpResidual()
{
  Real sum=0.0;
//  RealVectorType p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  for(unsigned int k=0; k<3; ++k)
    for(unsigned int l=0; l<3; ++l){
      sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(_component,_grad_test[_i][_qp],k,l)*p(k)*p(l);
    }
  return -0.5 * std::pow(_len_scale,2.0)*sum;
}

Real
FerroelectricCouplingU::computeQpJacobian()
{
  return 0.0;
}

Real
FerroelectricCouplingU::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component;
  Real sum=0.0;
//  RealVectorType w(_phi[_j][_qp]*_polar_x[_qp],_phi[_j][_qp]*_polar_y[_qp],_phi[_j][_qp]*_polar_z[_qp]);
  RealVectorValue w(_phi[_j][_qp] * _polar_x[_qp],_phi[_j][_qp] * _polar_y[_qp],_phi[_j][_qp] * _polar_z[_qp]);
  if( jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var){
    {
//    switch(jvar){
//    case _polar_x_var:
//      coupled_component=0;
//      break;
//    case _polar_y_var:
//      coupled_component=1;
//      break;
//    case _polar_z_var:
//      coupled_component=2;
//      break;
    if (jvar == _polar_x_var)
     {coupled_component=0;}
    else if (jvar == _polar_y_var)
     {coupled_component=1;}
    else if (jvar == _polar_z_var)
     {coupled_component=2;}
//    default:
    else
      mooseError("Something wrong with FerroelectricCoupling");
    }
    w(coupled_component)=w(coupled_component)*2.0;
    sum = _electrostrictive_tensor[_qp] . electrostrictiveProduct(_component, _grad_test[_i][_qp], coupled_component, w);
  }
  return  -0.5 * std::pow(_len_scale, 2.0) * sum;
}
