/**
 * @file   FerroelectricCouplingU.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15. 2015
 * @brief   Implement the kernel for displacement variable corresponding to ferroelectic coupling energy,
 *           Assume the energy has the form - 0.5 *q_ijkl * \frac{\partial u_i}{\partial x_j} * Pk_l where u is the displacement and P is the polarization.
 */

#include "FerroelectricCouplingU.h"

class FerroelectricCouplingU;

template<>
InputParameters validParams<FerroelectricCouplingU>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("disp_x", "The x component of the elasticity displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elasticity displacement");
  params.addRequiredCoupledVar("disp_z", "The z component of the elasticity displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

FerroelectricCouplingU::FerroelectricCouplingU(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _electrostrictive_tensor(getMaterialProperty<ElectrostrictiveTensorR4>("electrostrictive_tensor")),
   _component(getParam<unsigned int>("component")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FerroelectricCouplingU::computeQpResidual()
{
  Real sum = 0.0;
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(0, _polar_x_grad[_qp], _component, p);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(1, _polar_y_grad[_qp], _component, p);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(2, _polar_z_grad[_qp], _component, p);

  return std::pow(_len_scale, 3.0) * _test[_i][_qp] * sum;
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
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  if (jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
  {
    if (jvar == _polar_x_var)
    {
      coupled_component = 0;
    }
    else if (jvar == _polar_y_var)
    {
      coupled_component = 1;
    }
    else if (jvar == _polar_z_var)
    {
      coupled_component = 2;
    }
    sum1 = _electrostrictive_tensor[_qp].electrostrictiveProduct(_component, _grad_phi[_j][_qp], coupled_component, p);

    sum2 += _electrostrictive_tensor[_qp].electrostrictiveProduct(_component, _polar_x_grad[_qp], 0, coupled_component);
    sum2 += _electrostrictive_tensor[_qp].electrostrictiveProduct(_component, _polar_y_grad[_qp], 1, coupled_component);
    sum2 += _electrostrictive_tensor[_qp].electrostrictiveProduct(_component, _polar_z_grad[_qp], 2, coupled_component);

    return std::pow(_len_scale, 3.0) * (sum1 + sum2 * _phi[_j][_qp]) * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}
