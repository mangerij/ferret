/**
 * @file   PolarElectricPStrong.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @brief  PolarElectric interaction term;
 *
 */

#include "PolarElectricPStrong.h"

class PolarElectricPStrong;

template<>
InputParameters validParams<PolarElectricPStrong>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addCoupledVar("potential_ext", 0.0, "The external electric potential variable");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}



//Constructor
PolarElectricPStrong::PolarElectricPStrong(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_int_var(coupled("potential_int")),
   _potential_ext_var(coupled("potential_ext")),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
PolarElectricPStrong::computeQpResidual()
{
    Real RpolarP = 0.0;
    // TODO: Investigate the nature of this 1/2. This probably shouldn't be here...
    RpolarP += (0.5 * _potential_int_grad[_qp](_component) + _potential_ext_grad[_qp](_component)) * _test[_i][_qp] * std::pow(_len_scale, 2.0);

    //  Moose::out << "\n R_polarP-"; std::cout << _component << " = " << RpolarP;

    return RpolarP;
}

Real
PolarElectricPStrong::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricPStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
    if( jvar == _potential_int_var )
      return  0.5 * _grad_phi[_j][_qp](_component) * _test[_i][_qp] * std::pow(_len_scale, 2.0);
    else if( jvar == _potential_ext_var)
      return  _grad_phi[_j][_qp](_component) * _test[_i][_qp] * std::pow(_len_scale, 2.0);
    else
    {
      return 0.0;
    }
}
