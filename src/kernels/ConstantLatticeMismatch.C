/**
 * @file   ConstantLatticeMismatch.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 */

#include "ConstantLatticeMismatch.h"

template<>
InputParameters validParams<ConstantLatticeMismatch>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("Qxx", "dummy tensor xx");
  params.addRequiredCoupledVar("Qxy", "dummy tensor xy");
  params.addRequiredCoupledVar("Qxz", "dummy tensor xz");
  params.addRequiredCoupledVar("Qyx", "dummy tensor yx");
  params.addRequiredCoupledVar("Qyy", "dummy tensor yy");
  params.addRequiredCoupledVar("Qyz", "dummy tensor yz");
  params.addRequiredCoupledVar("Qzx", "dummy tensor zx");
  params.addRequiredCoupledVar("Qzy", "dummy tensor zy");
  params.addRequiredCoupledVar("Qzz", "dummy tensor zz");
  params.addRequiredCoupledVar("disp_x", "The x component of the displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the displacement");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<unsigned int>("deriv_component", "An integer corresponding to the direction of the derivative this kernel acts in. (0 for x, 1 for y, 2 for z)");
  return params;
}

ConstantLatticeMismatch::ConstantLatticeMismatch(const InputParameters & parameters)
  :Kernel(parameters),
   _Qxx(coupledValue("Qxx")),
   _Qxy(coupledValue("Qxy")),
   _Qxz(coupledValue("Qxz")),
   _Qyx(coupledValue("Qyx")),
   _Qyy(coupledValue("Qyy")),
   _Qyz(coupledValue("Qyz")),
   _Qzx(coupledValue("Qzx")),
   _Qzy(coupledValue("Qzy")),
   _Qzz(coupledValue("Qzz")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_y")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _component(getParam<unsigned int>("component")),
   _deriv_component(getParam<unsigned int>("deriv_component"))
{
}

Real
ConstantLatticeMismatch::computeQpResidual()
{
  if (_component == 0)
    {
      RealVectorValue Qx(_Qxx[_qp], _Qxy[_qp], _Qxz[_qp]);
      return _test[_i][_qp] * ( Qx(_deriv_component) - _disp_x_grad[_qp](_deriv_component));
    }
  else if (_component == 1)
    {
      RealVectorValue Qy(_Qyx[_qp], _Qyy[_qp], _Qyz[_qp]);
      return _test[_i][_qp] * ( Qy(_deriv_component) - _disp_y_grad[_qp](_deriv_component));
    }
  else if (_component == 2)
    {
      RealVectorValue Qz(_Qzx[_qp], _Qzy[_qp], _Qzz[_qp]);
      return _test[_i][_qp] * (Qz(_deriv_component) - _disp_z_grad[_qp](_deriv_component));
    }
  else
    {
      return 0.0;
    }
}

Real
ConstantLatticeMismatch::computeQpJacobian()
{
  if (_component == 0)
    {
      return _test[_i][_qp] * _phi[_j][_qp];
    }
  else if (_component == 1)
    {
      return _test[_i][_qp] * _phi[_j][_qp];
    }
  else if (_component == 2)
    {
      return _test[_i][_qp] * _phi[_j][_qp];
    }
  else
    {
      return 0.0;
    }
}

Real
ConstantLatticeMismatch::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _disp_x_var || jvar == _disp_y_var || jvar == _disp_z_var)
  {
    return - _test[_i][_qp] * _grad_phi[_j][_qp](_deriv_component);
  }
  else
  {
    return 0.0;
  }
}
