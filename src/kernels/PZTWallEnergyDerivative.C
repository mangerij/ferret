/**
 * @file   PZTWallEnergyDerivative.C
 * @author J. Manger <john.mangeri@uconn.edu>
 * @date   Dec 12 11:59:56 2015
 *
 * @brief
 *
 *
 */

#include "PZTWallEnergyDerivative.h"

template<>
InputParameters validParams<PZTWallEnergyDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("G110","Domain wall coefficient");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

PZTWallEnergyDerivative::PZTWallEnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
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
  _G110(getParam<Real>("G110")),
  _len_scale(getParam<Real>("len_scale"))
{
  //only for debug purpose
  std::cout<<"_G110="<<_G110<<"\n";
}

Real
PZTWallEnergyDerivative::computeQpResidual()
{
  Real Rwall = 0.0;

  Rwall += (_G110 * _polar_i_grad[_qp](_ii) * _grad_test[_i][_qp](_ii)) * _len_scale;

  //  Moose::out << "\n R_wall-"; std::cout << _component << " = " << Rwall;
  return Rwall;

}

Real
PZTWallEnergyDerivative::computeQpJacobian()
{
  return (_G110 * _grad_phi[_j][_qp](_ii) * _grad_test[_i][_qp](_ii) ) * _len_scale;
}

Real
PZTWallEnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  mooseAssert(jvar!=variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var)
  {
    const unsigned int _jj = (jvar==_polar_x_var)? 0: (jvar==_polar_y_var)? 1 : 2;
    return 0.0;
  }
  else
  {
    return 0.0;
  }
}
