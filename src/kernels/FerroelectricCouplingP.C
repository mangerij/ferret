/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "FerroelectricCouplingP.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ElectrostrictiveTensorTools.h"
#include "ComputeEigenstrain.h"
#include "libmesh/utility.h"

class FerroelectricCouplingP;

registerMooseObject("FerretApp", FerroelectricCouplingP);

template<>
InputParameters validParams<FerroelectricCouplingP>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t polarization of the electrostrictive coupling energy");
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elastic displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

FerroelectricCouplingP::FerroelectricCouplingP(const InputParameters & parameters)
  :Kernel(parameters),
   _electrostrictive_tensor(getMaterialProperty<RankFourTensor>("electrostrictive_tensor")),
   _eigenstrain(getMaterialProperty<RankTwoTensor>("eigenstrain")),
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
  Real sum = 0.0;
  Real RpCoupled = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  ///form three vectors of the _stress_free_strain[_qp] object

  RealVectorValue v0(_eigenstrain[_qp](0,0), _eigenstrain[_qp](0,1), _eigenstrain[_qp](0,2));
  RealVectorValue v1(_eigenstrain[_qp](1,0), _eigenstrain[_qp](1,1), _eigenstrain[_qp](1,2));
  RealVectorValue v2(_eigenstrain[_qp](2,0), _eigenstrain[_qp](2,1), _eigenstrain[_qp](2,2));
  

  //Moose::out << "\n e"; std::cout << 0; std::cout << 0; Moose::out << " = "; std::cout << _eigenstrain[_qp](0,0);
  //Moose::out << "\n e"; std::cout << 0; std::cout << 1; Moose::out << " = "; std::cout << _eigenstrain[_qp](0,1);
  //Moose::out << "\n e"; std::cout << 0; std::cout << 2; Moose::out << " = "; std::cout << _eigenstrain[_qp](0,2);
  //Moose::out << "\n e"; std::cout << 1; std::cout << 0; Moose::out << " = "; std::cout << _eigenstrain[_qp](1,0);
  //Moose::out << "\n e"; std::cout << 1; std::cout << 1; Moose::out << " = "; std::cout << _eigenstrain[_qp](1,1);
  //Moose::out << "\n e"; std::cout << 1; std::cout << 2; Moose::out << " = "; std::cout << _eigenstrain[_qp](1,2);
  //Moose::out << "\n e"; std::cout << 2; std::cout << 0; Moose::out << " = "; std::cout << _eigenstrain[_qp](2,0);
  //Moose::out << "\n e"; std::cout << 2; std::cout << 1; Moose::out << " = "; std::cout << _eigenstrain[_qp](2,1);
  //Moose::out << "\n e"; std::cout << 2; std::cout << 2; Moose::out << " = "; std::cout << _eigenstrain[_qp](2,2);

  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], _component, w);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], _component, w);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], _component, w);

  RpCoupled += Utility::pow<3>(_len_scale) * _test[_i][_qp] * sum; //factor of 1/2 disappears because variation of P_i P_j w.r.t P_k gives 2 terms equivalent in symmetry
  return - RpCoupled;
}

Real
FerroelectricCouplingP::computeQpJacobian()
{
  Real sum = 0.0;
  ///form three vectors of the _stress_free_strain[_qp] object
  RealVectorValue v0(_eigenstrain[_qp](0,0), _eigenstrain[_qp](0,1), _eigenstrain[_qp](0,2));
  RealVectorValue v1(_eigenstrain[_qp](1,0), _eigenstrain[_qp](1,1), _eigenstrain[_qp](1,2));
  RealVectorValue v2(_eigenstrain[_qp](2,0), _eigenstrain[_qp](2,1), _eigenstrain[_qp](2,2));
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], _component, _component);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], _component, _component);
  sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], _component, _component);
  return - Utility::pow<3>(_len_scale) * sum * _phi[_j][_qp] * _test[_i][_qp];
}

Real
FerroelectricCouplingP::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int coupled_component;
  Real sum = 0.0;
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  w(_component) = w(_component);

  ///form three vectors of the _stress_free_strain[_qp] object
  RealVectorValue v0(_eigenstrain[_qp](0,0), _eigenstrain[_qp](0,1), _eigenstrain[_qp](0,2));
  RealVectorValue v1(_eigenstrain[_qp](1,0), _eigenstrain[_qp](1,1), _eigenstrain[_qp](1,2));
  RealVectorValue v2(_eigenstrain[_qp](2,0), _eigenstrain[_qp](2,1), _eigenstrain[_qp](2,2));
  if( jvar == _polar_x_var || jvar == _polar_y_var || jvar == _polar_z_var)
  {
    if (jvar == _polar_x_var)
      {
        coupled_component = 0;

        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], _component, coupled_component);
      }
    else if (jvar == _polar_y_var)
      {
        coupled_component = 1;

        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], _component, coupled_component);
      }
    else if (jvar == _polar_z_var)
      {
        coupled_component = 2;

        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 0, v0 + _disp_x_grad[_qp], _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 1, v1 + _disp_y_grad[_qp], _component, coupled_component);
        sum += ElectrostrictiveTensorTools::electrostrictiveProduct(_electrostrictive_tensor[_qp], 2, v2 + _disp_z_grad[_qp], _component, coupled_component);
      }
    return - Utility::pow<3>(_len_scale) * sum * _phi[_j][_qp] * _test[_i][_qp];

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
    return - Utility::pow<3>(_len_scale) * sum * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}
