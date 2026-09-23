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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "CubicParentElasticADerivative.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "Assembly.h"

class CubicParentElasticADerivative;

registerMooseObject("FerretApp", CubicParentElasticADerivative);

InputParameters CubicParentElasticADerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t antiphase tilt of the rotostrictive coupling energy. Note: for cubic parent phase only.");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase (AFD) tilt");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase (AFD) tilt");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase (AFD) tilt");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addCoupledVar("displacements",
                       "The displacement variables. Optional; when given, the P<-u Jacobian block "
                       "is assembled (needs the 'elasticity_tensor' material property).");
  params.addCoupledVar("global_strain",
                       "The SCALAR mean strain, if it is being solved for. Optional; when given, the "
                       "P<->eps_bar blocks are assembled.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  return params;
}

CubicParentElasticADerivative::CubicParentElasticADerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _R11(getMaterialProperty<Real>("R11")),
   _R12(getMaterialProperty<Real>("R12")),
   _R44(getMaterialProperty<Real>("R44")),
   _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain")),
   _ndisp(coupledComponents("displacements")),
   _global_strain_var(isCoupledScalar("global_strain") ? coupledScalar("global_strain")
                                                       : libMesh::invalid_uint),
   _elasticity_tensor(_ndisp || _global_strain_var != libMesh::invalid_uint
                          ? &getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")
                          : nullptr)
{
  _disp_var.resize(_ndisp);
  for (unsigned int m = 0; m < _ndisp; ++m)
    _disp_var[m] = coupled("displacements", m);

}

Real
CubicParentElasticADerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return 2*(2*_C44[_qp]*_R44[_qp]*(_antiphase_A_x[_qp]*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_R44[_qp] - _antiphase_A_y[_qp]*_strain[_qp](0,1) - _antiphase_A_z[_qp]*_strain[_qp](0,2)) + 
     _C11[_qp]*_antiphase_A_x[_qp]*(2*Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp]*_R12[_qp] + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*Utility::pow<2>(_R12[_qp])) - _R11[_qp]*_strain[_qp](0,0) - 
        _R12[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](2,2)) + _C12[_qp]*_antiphase_A_x[_qp]*(4*Utility::pow<2>(_antiphase_A_x[_qp])*_R11[_qp]*_R12[_qp] + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + 
        Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) - 2*_R12[_qp]*_strain[_qp](0,0) - _R11[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](1,1) - _R11[_qp]*_strain[_qp](2,2) - _R12[_qp]*_strain[_qp](2,2)))*_test[_i][_qp];
  }
  else if (_component == 1)
  {
    return 2*(2*_C44[_qp]*_R44[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_R44[_qp] - _antiphase_A_x[_qp]*_strain[_qp](0,1) + _antiphase_A_z[_qp]*(_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_R44[_qp] - _strain[_qp](1,2))) + 
     _C11[_qp]*_antiphase_A_y[_qp]*(2*Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp]*_R12[_qp] + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*Utility::pow<2>(_R12[_qp])) - _R12[_qp]*_strain[_qp](0,0) - 
        _R11[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](2,2)) + _C12[_qp]*_antiphase_A_y[_qp]*(4*Utility::pow<2>(_antiphase_A_y[_qp])*_R11[_qp]*_R12[_qp] + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + 
        Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) - _R11[_qp]*_strain[_qp](0,0) - _R12[_qp]*_strain[_qp](0,0) - 2*_R12[_qp]*_strain[_qp](1,1) - _R11[_qp]*_strain[_qp](2,2) - _R12[_qp]*_strain[_qp](2,2)))*_test[_i][_qp];
  }
  else if (_component == 2)
  {
    return 2*(2*_C44[_qp]*_R44[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_R44[_qp] - _antiphase_A_x[_qp]*_strain[_qp](0,2) + _antiphase_A_y[_qp]*(_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_R44[_qp] - _strain[_qp](1,2))) + 
     _C11[_qp]*_antiphase_A_z[_qp]*(2*Utility::pow<2>(_antiphase_A_y[_qp])*_R11[_qp]*_R12[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*Utility::pow<2>(_R12[_qp])) - _R12[_qp]*_strain[_qp](0,0) - 
        _R12[_qp]*_strain[_qp](1,1) - _R11[_qp]*_strain[_qp](2,2)) + _C12[_qp]*_antiphase_A_z[_qp]*(4*Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp]*_R12[_qp] + 2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) - _R11[_qp]*_strain[_qp](0,0) - _R12[_qp]*_strain[_qp](0,0) - _R11[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](1,1) - 2*_R12[_qp]*_strain[_qp](2,2)))*_test[_i][_qp];
  }
  else
    return 0.0;
}

Real
CubicParentElasticADerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return 2*_phi[_j][_qp]*(2*_C44[_qp]*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_R44[_qp]) + _C11[_qp]*(2*Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp]*_R12[_qp] + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + 
        3*Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*Utility::pow<2>(_R12[_qp])) - _R11[_qp]*_strain[_qp](0,0) - _R12[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](2,2)) + 
     _C12[_qp]*(12*Utility::pow<2>(_antiphase_A_x[_qp])*_R11[_qp]*_R12[_qp] + 6*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + 
        Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) - 2*_R12[_qp]*_strain[_qp](0,0) - _R11[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](1,1) - _R11[_qp]*_strain[_qp](2,2) - _R12[_qp]*_strain[_qp](2,2)))*_test[_i][_qp];
  }
  else if (_component == 1)
  {
    return 2*_phi[_j][_qp]*(2*_C44[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_R44[_qp]) + _C11[_qp]*(2*Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp]*_R12[_qp] + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + 
        3*Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*Utility::pow<2>(_R12[_qp])) - _R12[_qp]*_strain[_qp](0,0) - _R11[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](2,2)) + 
     _C12[_qp]*(12*Utility::pow<2>(_antiphase_A_y[_qp])*_R11[_qp]*_R12[_qp] + 6*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + 
        Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) - _R11[_qp]*_strain[_qp](0,0) - _R12[_qp]*_strain[_qp](0,0) - 2*_R12[_qp]*_strain[_qp](1,1) - _R11[_qp]*_strain[_qp](2,2) - _R12[_qp]*_strain[_qp](2,2)))*_test[_i][_qp];
  }
  else if (_component == 2)
  {
    return 2*_phi[_j][_qp]*(2*_C44[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_R44[_qp]) + _C11[_qp]*(2*Utility::pow<2>(_antiphase_A_y[_qp])*_R11[_qp]*_R12[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + 
        3*Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*Utility::pow<2>(_R12[_qp])) - _R12[_qp]*_strain[_qp](0,0) - _R12[_qp]*_strain[_qp](1,1) - _R11[_qp]*_strain[_qp](2,2)) + 
     _C12[_qp]*(12*Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp]*_R12[_qp] + 6*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_R12[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + 
        Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) - _R11[_qp]*_strain[_qp](0,0) - _R12[_qp]*_strain[_qp](0,0) - _R11[_qp]*_strain[_qp](1,1) - _R12[_qp]*_strain[_qp](1,1) - 2*_R12[_qp]*_strain[_qp](2,2)))*_test[_i][_qp];
  }
  else
    return 0.0;
}

Real
CubicParentElasticADerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  // --- P <- u: only the strain in sigma = C:(eps - eps0) depends on u ---
  for (unsigned int m = 0; m < _ndisp; ++m)
    if (jvar == _disp_var[m])
    {
      RankTwoTensor deps;   // sym(e_m (x) grad phi)
      for (unsigned int d = 0; d < 3; ++d)
      {
        deps(m, d) += 0.5 * _grad_phi[_j][_qp](d);
        deps(d, m) += 0.5 * _grad_phi[_j][_qp](d);
      }
      const RankTwoTensor dsigma = (*_elasticity_tensor)[_qp] * deps;
      return -dsigma.doubleContraction(dEigenstrain_dA(_component)) * _test[_i][_qp];
    }

  if (_component == 0)
  {
    if (jvar == _antiphase_A_y_var)
    {
      return 4*_phi[_j][_qp]*(_C11[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + _C12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + _C44[_qp]*_R44[_qp]*(2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_R44[_qp] - _strain[_qp](0,1)))*_test[_i][_qp];
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return 4*_phi[_j][_qp]*(_C11[_qp]*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + _C12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + _C44[_qp]*_R44[_qp]*(2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_R44[_qp] - _strain[_qp](0,2)))*_test[_i][_qp];
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiphase_A_x_var)
    {
      return 4*_phi[_j][_qp]*(_C11[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + _C12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + _C44[_qp]*_R44[_qp]*(2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_R44[_qp] - _strain[_qp](0,1)))*_test[_i][_qp];
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return 4*_phi[_j][_qp]*(_C11[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + _C12[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + _C44[_qp]*_R44[_qp]*(2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_R44[_qp] - _strain[_qp](1,2)))*_test[_i][_qp];
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiphase_A_x_var)
    {
      return 4*_phi[_j][_qp]*(_C11[_qp]*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + _C12[_qp]*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + _C44[_qp]*_R44[_qp]*(2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_R44[_qp] - _strain[_qp](0,2)))*_test[_i][_qp];
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return 4*_phi[_j][_qp]*(_C11[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_R12[_qp]*(2*_R11[_qp] + _R12[_qp]) + _C12[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*(Utility::pow<2>(_R11[_qp]) + 2*_R11[_qp]*_R12[_qp] + 3*Utility::pow<2>(_R12[_qp])) + _C44[_qp]*_R44[_qp]*(2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_R44[_qp] - _strain[_qp](1,2)))*_test[_i][_qp];
    }
    else
    {
      return 0.0;
    }
  }
  return 0.0;
}


RankTwoTensor
CubicParentElasticADerivative::dEigenstrain_dA(unsigned int i) const
{
  const Real P[3] = {_antiphase_A_x[_qp], _antiphase_A_y[_qp], _antiphase_A_z[_qp]};
  RankTwoTensor D;
  for (unsigned int j = 0; j < 3; ++j)
  {
    if (j == i)
      D(j, j) = 2.0 * _R11[_qp] * P[i];
    else
    {
      D(j, j) = 2.0 * _R12[_qp] * P[i];
      D(i, j) = D(j, i) = _R44[_qp] * P[j];
    }
  }
  return D;
}

RankTwoTensor
CubicParentElasticADerivative::unitStrain(unsigned int b)
{
  RankTwoTensor E;
  switch (b)
  {
    case 0: E(0, 0) = 1; break;
    case 1: E(1, 1) = 1; break;
    case 2: E(2, 2) = 1; break;
    case 3: E(1, 2) = E(2, 1) = 1; break;
    case 4: E(0, 2) = E(2, 0) = 1; break;
    case 5: E(0, 1) = E(1, 0) = 1; break;
  }
  return E;
}

void
CubicParentElasticADerivative::computeOffDiagJacobianScalar(unsigned int jvar)
{
  // P <-> eps_bar. The mean strain enters sigma exactly like a uniform strain, so
  //     dR_i / d eps_bar_b = -(C:E_b) : d eps0/dP_i  psi ,
  // and the eps_bar residual (integral of sigma, MOOSE's GlobalStrainUserObject) has
  //     d R_eps_bar_b / dP_i = -(C:E_b) : d eps0/dP_i  phi ,
  // the exact transpose. MOOSE's GlobalStrain scalar kernel supplies only its own diagonal,
  // so both blocks are assembled here; each P component fills its own row, nothing twice.
  if (_global_strain_var == libMesh::invalid_uint || jvar != _global_strain_var)
    return;

  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);

  prepareMatrixTag(_assembly, _var.number(), jvar);
  for (_i = 0; _i < _test.size(); ++_i)
    for (_j = 0; _j < jv.order(); ++_j)
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      {
        const RankTwoTensor dsigma = (*_elasticity_tensor)[_qp] * unitStrain(_j);
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] *
                             (-dsigma.doubleContraction(dEigenstrain_dA(_component))) *
                             _test[_i][_qp];
      }
  accumulateTaggedLocalMatrix();

  // Transpose for the eps_bar rows. MOOSE's scalar residual is int sigma_ij dV for the SINGLE
  // tensor component (ij) of each scalar entry, whereas the unit strain E_b used above moves
  // both eps_ij and eps_ji: the column block is right as it is, the row block for the three
  // shear components is exactly half of the transpose.
  _ke_copy = _local_ke;
  prepareMatrixTag(_assembly, jvar, _var.number());
  _ke_copy.get_transpose(_local_ke);
  for (unsigned int b = 3; b < _local_ke.m(); ++b)
    for (unsigned int c = 0; c < _local_ke.n(); ++c)
      _local_ke(b, c) *= 0.5;
  accumulateTaggedLocalMatrix();
}
