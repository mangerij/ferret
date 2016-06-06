/**
 * @file   FerroelectricCouplingP.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15. 2015
 *
 * @brief  Implement the kernel for polar variables corresponding to ferroelectic coupling energy after
 *         the variational derivative of the polar dependent terms have been taken.
 *         This is only the quadratic term. See notes.
 */

#include "FerroelectricCouplingP.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ElectrostrictiveTensorTools.h"
#include "ComputeEigenstrain.h"

class FerroelectricCouplingP;

template<>
InputParameters validParams<FerroelectricCouplingP>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0,  "The z component of the elastic displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("strain_scale", 1.0, "the strain_scale");
  params.addParam<Real>("artificial", 1.0, "artificial increase coupling");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

FerroelectricCouplingP::FerroelectricCouplingP(const InputParameters & parameters)
  :Kernel(parameters),
   _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
   _stress_free_strain(getMaterialProperty<RankTwoTensor>("stress_free_strain")),
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
   _strain_scale(getParam<Real>("strain_scale")),
   _artificial(getParam<Real>("artificial")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FerroelectricCouplingP::computeQpResidual()
{
  Real sum = 0.0;
  Real RpCoupled = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  //form three vectors of the _stress_free_strain[_qp] object

  RealVectorValue v0(_stress_free_strain[_qp](0,0), _stress_free_strain[_qp](0,1), _stress_free_strain[_qp](0,2));
  RealVectorValue v1(_stress_free_strain[_qp](1,0), _stress_free_strain[_qp](1,1), _stress_free_strain[_qp](1,2));
  RealVectorValue v2(_stress_free_strain[_qp](2,0), _stress_free_strain[_qp](2,1), _stress_free_strain[_qp](2,2));

  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, _component, w);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, _component, w);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, _component, w);

  RpCoupled += std::pow(_len_scale, 3.0) * _test[_i][_qp] * sum;
  return - _artificial * RpCoupled;
}

Real
FerroelectricCouplingP::computeQpJacobian()
{
  Real sum = 0.0;
  //form three vectors of the _stress_free_strain[_qp] object
  RealVectorValue v0(_stress_free_strain[_qp](0,0), _stress_free_strain[_qp](0,1), _stress_free_strain[_qp](0,2));
  RealVectorValue v1(_stress_free_strain[_qp](1,0), _stress_free_strain[_qp](1,1), _stress_free_strain[_qp](1,2));
  RealVectorValue v2(_stress_free_strain[_qp](2,0), _stress_free_strain[_qp](2,1), _stress_free_strain[_qp](2,2));
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, _component, _component);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, _component, _component);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, _component, _component);
  return - _artificial * std::pow(_len_scale, 3.0) * sum * _phi[_j][_qp] * _test[_i][_qp];
}

Real
FerroelectricCouplingP::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component;
  Real sum = 0.0;
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  w(_component) = w(_component);

  //form three vectors of the _stress_free_strain[_qp] object
  RealVectorValue v0(_stress_free_strain[_qp](0,0), _stress_free_strain[_qp](0,1), _stress_free_strain[_qp](0,2));
  RealVectorValue v1(_stress_free_strain[_qp](1,0), _stress_free_strain[_qp](1,1), _stress_free_strain[_qp](1,2));
  RealVectorValue v2(_stress_free_strain[_qp](2,0), _stress_free_strain[_qp](2,1), _stress_free_strain[_qp](2,2));
  if( jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
  {
    if (jvar == _polar_x_var)
      {
        coupled_component = 0;

        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, _component, coupled_component);
      }
    else if (jvar == _polar_y_var)
      {
        coupled_component = 1;

        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, _component, coupled_component);
      }
    else if (jvar == _polar_z_var)
      {
        coupled_component = 2;

        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, _disp_x_grad[_qp] - v0, _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, _disp_y_grad[_qp] - v1, _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, _disp_z_grad[_qp] - v2, _component, coupled_component);
      }
    return - _artificial * std::pow(_len_scale, 3.0) * sum * _phi[_j][_qp] * _test[_i][_qp];

  }
  else if(jvar == _disp_x_var || jvar == _disp_y_var || jvar == _disp_z_var)
  {
    if (jvar == _disp_x_var)
      {
        coupled_component = 0;
        sum = ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], coupled_component, _grad_phi[_j][_qp], _component, w);
      }
    else if (jvar == _disp_y_var)
      {
        coupled_component = 1;
        sum = ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], coupled_component, _grad_phi[_j][_qp], _component, w);
      }
    else if (jvar == _disp_z_var)
      {
        coupled_component = 2;
        sum = ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], coupled_component, _grad_phi[_j][_qp], _component, w);
      }
    return - _artificial * std::pow(_len_scale, 3.0) * sum * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}
