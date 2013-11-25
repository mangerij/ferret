/**
 * @file   PolarElectricEStrong.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 11:59:56 2013
 *
 * @brief
 *
 *
 */

#include "PolarElectricEStrong.h"

class PolarElectricEStrong;

template<>
InputParameters validParams<PolarElectricEStrong>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  params.addParam<Real>("polar_electric_scale",1.0,"polar_electric scale");
  return params;
}



//Constructor
PolarElectricEStrong::PolarElectricEStrong(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _permittivity(getParam<Real>("permittivity")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _polar_electric_scale(getParam<Real>("polar_electric_scale"))
{}



Real
PolarElectricEStrong::computeQpResidual()
{
  return -(_polar_x[_qp]*_grad_test[_i][_qp](0)+_polar_y[_qp]*_grad_test[_i][_qp](1)+_polar_z[_qp]*_grad_test[_i][_qp](2))*pow(_len_scale,2.0)*_polar_electric_scale;
}

Real
PolarElectricEStrong::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricEStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
  if( jvar == coupled("polar_x") )
    return -_phi[_j][_qp]*_grad_test[_i][_qp](0)*pow(_len_scale,2.0)*_polar_electric_scale;
  else if( jvar == coupled("polar_y"))
    return -_phi[_j][_qp]*_grad_test[_i][_qp](1)*pow(_len_scale,2.0)*_polar_electric_scale;
  else if(jvar == coupled("polar_z"))
    return -_phi[_j][_qp]*_grad_test[_i][_qp](2)*pow(_len_scale,2.0)*_polar_electric_scale;
  else{
    return 0.0;
  }
}
