/**
 * @file   FerroelectricCouplingP.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15. 2015
 *
 * @brief  Implement the kernel for polar variables corresponding to ferroelectic coupling energy.
 *         Assume the energy has the form -0.5*q_ijkl* ui_j * Pk*Pl where u is the displacement and P is the polarization.
 */

#include "FerroelectricCouplingP.h"

class FerroelectricCouplingP;

template<>
InputParameters validParams<FerroelectricCouplingP>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "The x component of the elasticity displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elasticity displacement");
  params.addRequiredCoupledVar("disp_z", "The z component of the elasticity displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

FerroelectricCouplingP::FerroelectricCouplingP(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _electrostrictive_tensor(getMaterialProperty<ElectrostrictiveTensorR4>("electrostrictive_tensor")),
   _component(getParam<unsigned int>("component")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FerroelectricCouplingP::computeQpResidual()
{
  // Compute terms of -0.5 q_{ijkl} \epsilon_{ij} P_k \delta P_l for l = _component and
  // \epsilon_{ij} = 0.5 (\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i})
  // the minor symmetry of q_{ijkl = jikl} allows the below sum to only care about one disp_grad
  // vector and not its transpose as well.

  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  Real sum = 0.0;
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(0, _disp_x_grad[_qp], _component, w);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(1, _disp_y_grad[_qp], _component, w);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(2, _disp_z_grad[_qp], _component, w);
  return - 0.5 * std::pow(_len_scale, 3.0) * _test[_i][_qp] * sum;
}

Real
FerroelectricCouplingP::computeQpJacobian()
{
  Real sum = 0.0;
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(0, _disp_x_grad[_qp], _component, _component);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(1, _disp_y_grad[_qp], _component, _component);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(2, _disp_z_grad[_qp], _component, _component);
  return - 0.5 * std::pow(_len_scale, 3.0) * sum * _phi[_j][_qp] * _test[_i][_qp];
}

Real
FerroelectricCouplingP::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component;
  Real sum=0.0;

  if( jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
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
  }
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(0, _disp_x_grad[_qp], _component, coupled_component);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(1, _disp_y_grad[_qp], _component, coupled_component);
  sum += _electrostrictive_tensor[_qp].electrostrictiveProduct(2, _disp_z_grad[_qp], _component, coupled_component);
  return - 0.5 * std::pow(_len_scale, 3.0) * sum * _phi[_j][_qp] * _test[_i][_qp];
}
