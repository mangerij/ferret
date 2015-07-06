/**
 * @file   FerroelectricCouplingUWeak.C
 * @author J. Mangeri <mangerij@anl.gov>
 * @date   Jul. 1. 2015
 * @brief   Implement the kernel for displacement variable corresponding to ferroelectic coupling energy,
 *          Assume the energy has the form -0.5*q_ijkl* ui_j * Pk_l where u is the displacement and P is the polarization.
 *          with zero jacobian terms (weak coupling)
 */

#include "FerroelectricCouplingUWeak.h"

class FerroelectricCouplingUWeak;

template<>
InputParameters validParams<FerroelectricCouplingUWeak>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

FerroelectricCouplingUWeak::FerroelectricCouplingUWeak(const std::string & name, InputParameters parameters)
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
FerroelectricCouplingUWeak::computeQpResidual()
{
  Real sum=0.0;
//  RealVectorType p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  for(unsigned int k=0; k<3; ++k)
    for(unsigned int l=0; l<3; ++l)
    {
      sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(_component, _grad_test[_i][_qp], k, l) * p(k) * p(l);
    }
  return - 0.5 * std::pow(_len_scale, 2.0) * sum;
}

Real
FerroelectricCouplingUWeak::computeQpJacobian()
{
  return 0.0;
}

Real
FerroelectricCouplingUWeak::computeQpOffDiagJacobian(unsigned int jvar)
{
  return  0.0;
}
