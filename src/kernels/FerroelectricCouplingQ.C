/**
 * @file   FerroelectricCouplingQ.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15. 2015
 *
 * @brief  Implement the kernel for polar variables corresponding to ferroelectic coupling energy after
 *         the variational derivative of the polar dependent terms have been taken.
 *         This is only the quartic term. See notes.
 */

#include "FerroelectricCouplingQ.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ElectrostrictiveTensorTools.h"

class FerroelectricCouplingQ;

template<>
InputParameters validParams<FerroelectricCouplingQ>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("artificial", 1.0, "artificial increase coupling");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

FerroelectricCouplingQ::FerroelectricCouplingQ(const InputParameters & parameters)
  :Kernel(parameters),
   _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
   _electrostrictive_tensorQ(getMaterialProperty<RankFourTensor>("electrostrictive_tensorQ")),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_y")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _artificial(getParam<Real>("artificial")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FerroelectricCouplingQ::computeQpResidual()
{
  Real sum = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 0, w, _component, w) * w(0);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 1, w, _component, w) * w(1);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 2, w, _component, w) * w(2);
  return 2.0 * std::pow(_len_scale, 3.0) * _test[_i][_qp] * sum;
}

Real
FerroelectricCouplingQ::computeQpJacobian()
{
  Real sum = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 0, w, _component, _component) * w(0);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 1, w, _component, _component) * w(1);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 2, w, _component, _component) * w(2);
  return 2.0 * std::pow(_len_scale, 3.0) * _phi[_j][_qp] * _test[_i][_qp] * sum;
}

Real
FerroelectricCouplingQ::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component;
  Real sum = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  if( jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
  {
    if (jvar == _polar_x_var)
      {
        coupled_component = 0;
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 0, w, _component, coupled_component) * w(0);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 1, w, _component, coupled_component) * w(1);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 2, w, _component, coupled_component) * w(2);
      }
    else if (jvar == _polar_y_var)
      {
        coupled_component = 1;
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 0, w, _component, coupled_component) * w(0);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 1, w, _component, coupled_component) * w(1);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 2, w, _component, coupled_component) * w(2);
      }
    else if (jvar == _polar_z_var)
      {
        coupled_component = 2;
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 0, w, _component, coupled_component) * w(0);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 1, w, _component, coupled_component) * w(1);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensorQ[_qp], 2, w, _component, coupled_component) * w(2);
      }
    return 2.0 * std::pow(_len_scale, 3.0) * sum * _phi[_j][_qp] * _test[_i][_qp];

  }
  else
  {
    return 0.0;
  }
}
