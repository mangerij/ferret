/**
 * @file   PolarElectricPStrong.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 11:59:56 2013
 *
 * @brief  PolarElectric interaction term,weak coupling-- Polar part;
 *
 *
 */

#include "PolarElectricPStrong.h"

class PolarElectricPStrong;

template<>
InputParameters validParams<PolarElectricPStrong>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<Real>("permittivity_ext", "external permittivity");
  params.addRequiredParam<Real>("permittivity_int", "internal permittivity");
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  params.addParam<Real>("energy_scale",1.0,"energy scale");
  return params;
}



//Constructor
PolarElectricPStrong::PolarElectricPStrong(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _component(getParam<unsigned int>("component")),
   _permittivity_int(getParam<Real>("permittivity_int")),
   _permittivity_ext(getParam<Real>("permittivity_ext")),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext")),
   _len_scale(getParam<Real>("len_scale")),
   _energy_scale(getParam<Real>("energy_scale"))
{}



Real
PolarElectricPStrong::computeQpResidual()
{
    return (0.5*_permittivity_int*_potential_int_grad[_qp](_component)+_permittivity_ext*_potential_ext_grad[_qp](_component))*_test[_i][_qp]*pow(_len_scale,2.0)*_energy_scale;
}

Real
PolarElectricPStrong::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricPStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
  if( jvar == coupled("potential_int") )
    return _permittivity_int*0.5*_grad_phi[_j][_qp](_component)*_test[_i][_qp]*pow(_len_scale,2.0)*_energy_scale;
  else if( jvar == coupled("potential_ext"))
    return _permittivity_ext*_grad_phi[_j][_qp](_component)*_test[_i][_qp]*pow(_len_scale,2.0)*_energy_scale;
  else{
    return 0.0;
  }
}
