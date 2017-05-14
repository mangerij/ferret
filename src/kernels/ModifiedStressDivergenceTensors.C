/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "StressDivergenceTensors.h"
#include "ModifiedStressDivergenceTensors.h"
#include "Material.h"
#include "MooseMesh.h"
#include "ElasticityTensorTools.h"

template<>
InputParameters validParams<ModifiedStressDivergenceTensors>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Stress divergence kernel (used by the TensorMechanics action)");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements", "The string of displacements suitable for the problem statement");
  params.addCoupledVar("temp", "The temperature");
  params.addParam<std::string>("base_name", "Material property base name");
  params.set<bool>("use_displaced_mesh") = false;
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}


ModifiedStressDivergenceTensors::ModifiedStressDivergenceTensors(const InputParameters & parameters) :
    Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(3),
    _temp_coupled(isCoupled("temp")),
    _temp_var(_temp_coupled ? coupled("temp") : 0),
    _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
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
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension.");
}

Real
ModifiedStressDivergenceTensors::computeQpResidual()
{
  Real sum = 0.0;
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], 0, p) * _polar_x[_qp];
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], 1, p) * _polar_y[_qp];
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], 2, p) * _polar_z[_qp];
  return (_stress[_qp].row(_component)) * _grad_test[_i][_qp] - sum;
}

Real
ModifiedStressDivergenceTensors::computeQpJacobian()
{
  return ElasticityTensorTools::elasticJacobian(_Jacobian_mult[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);
}

Real
ModifiedStressDivergenceTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component = 0;
  bool active(false);
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  unsigned int coupled_component1;
  RealVectorValue p(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  // return = 0.0;
  for (unsigned int i = 0; i < _ndisp; ++i)
    if (jvar == _disp_var[i])
    {
      coupled_component = i;
      active = true;
    }

  if (active)
    return ElasticityTensorTools::elasticJacobian(_Jacobian_mult[_qp], _component, coupled_component,
                                          _grad_test[_i][_qp], _grad_phi[_j][_qp]) ;
  if (_temp_coupled && jvar == _temp_var)
  {
    //return _d_stress_dT[_qp].rowDot(_component, _grad_test[_i][_qp]) * _phi[_j][_qp];
    return 0.0;
  }
  if (jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
  {
    if (jvar == _polar_x_var)
    {
      coupled_component1 = 0;
      sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], coupled_component1, p);
    }
    else if (jvar == _polar_y_var)
    {
      coupled_component1 = 1;
      sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], coupled_component1, p);
    }
    else if (jvar == _polar_z_var)
    {
      coupled_component1 = 2;
      sum1 += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], _component, _grad_test[_i][_qp], coupled_component1, p) ;
    }
    return - 2.0 * std::pow(_len_scale, 2.0) * _phi[_j][_qp] * sum1;
  }
  else
  {
    return 0.0;
  }
}
