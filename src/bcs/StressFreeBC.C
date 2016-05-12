/************************************************************************/
/* StressFree BC:                                                       */
/*     This BC is intended to implement sigma_i = 0                     */
/*     at a boundary for i = 0,1,2                                      */
/************************************************************************/

#include "StressFreeBC.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ElectrostrictiveTensorTools.h"

template<>
InputParameters validParams<StressFreeBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
    params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
    params.addCoupledVar("disp_z", 0.0,  "The z component of the elastic displacement");
    params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
    params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
    params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
    params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

StressFreeBC::StressFreeBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
  _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
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
  _disp_x_grad(coupledGradient("disp_x")),
  _disp_y_grad(coupledGradient("disp_y")),
  _disp_z_grad(coupledGradient("disp_z"))
{
}

Real
StressFreeBC::computeQpResidual()
{
  Real sum1 = 0;
  Real sum2 = 0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, w, _component, w); //contracts w^2 on q_ijkl
  sum2 += _elasticity_tensor[_qp](_component, 2, 0, 0) * _disp_x_grad[_qp](0) + _elasticity_tensor[_qp](_component, 2, 1, 0) * _disp_y_grad[_qp](0)
          + _elasticity_tensor[_qp](_component, 2, 2, 0) * _disp_z_grad[_qp](0) + _elasticity_tensor[_qp](_component, 2, 0, 1) * _disp_x_grad[_qp](1)
          + _elasticity_tensor[_qp](_component, 2, 1, 1) * _disp_y_grad[_qp](1) + _elasticity_tensor[_qp](_component, 2, 2, 1) * _disp_z_grad[_qp](1)
          + _elasticity_tensor[_qp](_component, 2, 0, 2) * _disp_x_grad[_qp](2) + _elasticity_tensor[_qp](_component, 2, 1, 2) * _disp_y_grad[_qp](2)
          + _elasticity_tensor[_qp](_component, 2, 2, 2) * _disp_z_grad[_qp](2);

  return  (sum2 - sum1) *_test[_i][_qp];
}

Real
StressFreeBC::computeQpJacobian()
{
    Real sum3 = 0.0;
    if (_component == 0)
      {
        sum3 += _elasticity_tensor[_qp](_component, 2, 0, 0) * _grad_phi[_j][_qp](0) + _elasticity_tensor[_qp](_component, 2, 0, 1) * _grad_phi[_j][_qp](1)
                + _elasticity_tensor[_qp](_component, 2, 0, 2) * _grad_phi[_j][_qp](2);
      }
    else if (_component == 1)
      {
        sum3 += _elasticity_tensor[_qp](_component, 2, 1, 0) * _grad_phi[_j][_qp](0) + _elasticity_tensor[_qp](_component, 2, 1, 1) * _grad_phi[_j][_qp](1)
                + _elasticity_tensor[_qp](_component, 2, 1, 2) * _grad_phi[_j][_qp](2);
      }
    else if (_component == 2)
      {
        sum3 += _elasticity_tensor[_qp](_component, 2, 2, 0) * _grad_phi[_j][_qp](0) + _elasticity_tensor[_qp](_component, 2, 2, 1) * _grad_phi[_j][_qp](1)
                + _elasticity_tensor[_qp](_component, 2, 2, 2) * _grad_phi[_j][_qp](2);
      }
      return _test[_i][_qp] * sum3;
}

Real
StressFreeBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real sum1 = 0;
  Real sum2 = 0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  unsigned int coupled_component;
  if( jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
  {
    if (jvar == _polar_x_var)
      {
        coupled_component = 0;
        sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, w, _component, coupled_component);
      }
    else if (jvar == _polar_y_var)
      {
        coupled_component = 1;
        sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, w, _component, coupled_component);
      }
    else if (jvar == _polar_z_var)
      {
        coupled_component = 2;
        sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, w, _component, coupled_component);
      }
    return - sum1 * _phi[_j][_qp] * _test[_i][_qp];

  }
  else if(jvar == _disp_x_var || jvar == _disp_y_var || jvar == _disp_z_var)
  {
    if (jvar == _disp_x_var)
      {
        coupled_component = 0;
        sum2 += _elasticity_tensor[_qp](_component, 2, 0, 0) * _grad_phi[_j][_qp](0) + _elasticity_tensor[_qp](_component, 2, 0, 1) * _grad_phi[_j][_qp](1)
                + _elasticity_tensor[_qp](_component, 2, 0, 2) * _grad_phi[_j][_qp](2);
      }
    else if (jvar == _disp_y_var)
      {
        coupled_component = 1;
        sum2 += _elasticity_tensor[_qp](_component, 2, 1, 0) * _grad_phi[_j][_qp](0) + _elasticity_tensor[_qp](_component, 2, 1, 1) * _grad_phi[_j][_qp](1)
                + _elasticity_tensor[_qp](_component, 2, 1, 2) * _grad_phi[_j][_qp](2);
      }
    else if (jvar == _disp_z_var)
      {
        coupled_component = 2;
        sum2 += _elasticity_tensor[_qp](_component, 2, 2, 0) * _grad_phi[_j][_qp](0) + _elasticity_tensor[_qp](_component, 2, 2, 1) * _grad_phi[_j][_qp](1)
                + _elasticity_tensor[_qp](_component, 2, 2, 2) * _grad_phi[_j][_qp](2);
      }
    return sum2 * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}
