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

#include "CubicParentFlexoelectricFullPDerivative.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "Assembly.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CubicParentFlexoelectricFullPDerivative);

/*
  Generated from ONE energy density (scratch_calcs/flexo_port/gen.py),

     W = 1/2 (eps - eps^Q - eps^F) : C : (eps - eps^Q - eps^F),

  so the residual, every Jacobian block and the energy postprocessor are exact derivatives
  of each other.  This kernel carries BOTH the electrostrictive and the flexoelectric parts
  of dW/dP: use it INSTEAD of CubicParentElasticPDerivative, never together with it.

  The flexoelectric eigenstrain is eps^F_ij = -A_ijkl dP_k/dx_l with A = s : f, where the
  flexocoupling tensor f is entered through its TENSOR components mu11 = f_1111,
  mu12 = f_1122, mu44 = f_1212 and the compliances in Voigt form (s_1212 = s44/4).
  With this sign, the published f = C : F of Gu et al., PRB 89, 174111 (2014), entered
  directly, gives their tail-to-tail Neel component.  The eigenstrain enters the
  CompositeEigenstrain with weight +1.
*/

InputParameters CubicParentFlexoelectricFullPDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual of the full elastic energy, electrostrictive AND flexoelectric, w.r.t. the polarization, for a cubic parent phase. Replaces CubicParentElasticPDerivative when flexoelectricity is on.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addCoupledVar("displacements",
                       "The displacement variables. Optional; when given, the P<-u Jacobian block "
                       "is assembled.");
  params.addCoupledVar("global_strain",
                       "The SCALAR mean strain, if it is being solved for. Optional; when given, the "
                       "P<->eps_bar blocks are assembled.");
  return params;
}

CubicParentFlexoelectricFullPDerivative::CubicParentFlexoelectricFullPDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _Q11(getMaterialProperty<Real>("Q11")),
   _Q12(getMaterialProperty<Real>("Q12")),
   _Q44(getMaterialProperty<Real>("Q44")),
   _s11(getMaterialProperty<Real>("s11")),
   _s12(getMaterialProperty<Real>("s12")),
   _s44(getMaterialProperty<Real>("s44")),
   _mu11(getMaterialProperty<Real>("mu11")),
   _mu12(getMaterialProperty<Real>("mu12")),
   _mu44(getMaterialProperty<Real>("mu44")),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
   _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain")),
   _ndisp(coupledComponents("displacements")),
   _global_strain_var(isCoupledScalar("global_strain") ? coupledScalar("global_strain")
                                                       : libMesh::invalid_uint)
{
  _disp_var.resize(_ndisp);
  for (unsigned int m = 0; m < _ndisp; ++m)
    _disp_var[m] = coupled("displacements", m);
}

Real
CubicParentFlexoelectricFullPDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](1)*_polar_y[_qp]*_test[_i][_qp] + _polar_x_grad[_qp](2)*_polar_z[_qp]*_test[_i][_qp] + _polar_y_grad[_qp](0)*_polar_y[_qp]*_test[_i][_qp]
      + _polar_z_grad[_qp](0)*_polar_z[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](1)*_polar_x[_qp]*_polar_y[_qp] + _grad_test[_i][_qp](2)*_polar_x[_qp]*_polar_z[_qp])
      - 4*_C44[_qp]*_Q44[_qp]*(_strain[_qp](0,1)*_polar_y[_qp]*_test[_i][_qp] + _strain[_qp](0,2)*_polar_z[_qp]*_test[_i][_qp]) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(_polar_x_grad[_qp](1)*_grad_test[_i][_qp](1)
      + _polar_x_grad[_qp](2)*_grad_test[_i][_qp](2) + _polar_y_grad[_qp](0)*_grad_test[_i][_qp](1) + _polar_z_grad[_qp](0)*_grad_test[_i][_qp](2)) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(_strain[_qp](0,1)*_grad_test[_i][_qp](1)
      + _strain[_qp](0,2)*_grad_test[_i][_qp](2)) + _strain[_qp](0,0)*_grad_test[_i][_qp](0)*(_C11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_strain[_qp](0,0)*_polar_x[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]) + _polar_x_grad[_qp](0)*_grad_test[_i][_qp](0)*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*_polar_x_grad[_qp](0)*_polar_x[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) - _grad_test[_i][_qp](0)*Utility::pow<2>(_polar_x[_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*Utility::pow<3>(_polar_x[_qp])*_test[_i][_qp]*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + (_strain[_qp](1,1)*_grad_test[_i][_qp](0)
      + _strain[_qp](2,2)*_grad_test[_i][_qp](0))*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp]
      + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (_polar_y_grad[_qp](1)*_grad_test[_i][_qp](0) + _polar_z_grad[_qp](2)*_grad_test[_i][_qp](0))*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - (_grad_test[_i][_qp](0)*Utility::pow<2>(_polar_y[_qp]) + _grad_test[_i][_qp](0)*Utility::pow<2>(_polar_z[_qp]))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*(_strain[_qp](1,1)*_polar_x[_qp]*_test[_i][_qp] + _strain[_qp](2,2)*_polar_x[_qp]*_test[_i][_qp])*(_C11[_qp]*_Q12[_qp]
      + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_polar_y_grad[_qp](1)*_polar_x[_qp]*_test[_i][_qp] + _polar_z_grad[_qp](2)*_polar_x[_qp]*_test[_i][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_test[_i][_qp]
      + _polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_test[_i][_qp])*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp] + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp]) + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp]));
  }
  else if (_component == 1)
  {
    return -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](1)*_polar_x[_qp]*_test[_i][_qp] + _polar_y_grad[_qp](0)*_polar_x[_qp]*_test[_i][_qp] + _polar_y_grad[_qp](2)*_polar_z[_qp]*_test[_i][_qp]
      + _polar_z_grad[_qp](1)*_polar_z[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](0)*_polar_x[_qp]*_polar_y[_qp] + _grad_test[_i][_qp](2)*_polar_y[_qp]*_polar_z[_qp])
      - 4*_C44[_qp]*_Q44[_qp]*(_strain[_qp](0,1)*_polar_x[_qp]*_test[_i][_qp] + _strain[_qp](1,2)*_polar_z[_qp]*_test[_i][_qp]) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(_polar_x_grad[_qp](1)*_grad_test[_i][_qp](0)
      + _polar_y_grad[_qp](0)*_grad_test[_i][_qp](0) + _polar_y_grad[_qp](2)*_grad_test[_i][_qp](2) + _polar_z_grad[_qp](1)*_grad_test[_i][_qp](2)) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(_strain[_qp](0,1)*_grad_test[_i][_qp](0)
      + _strain[_qp](1,2)*_grad_test[_i][_qp](2)) + _strain[_qp](1,1)*_grad_test[_i][_qp](1)*(_C11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_strain[_qp](1,1)*_polar_y[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]) + _polar_y_grad[_qp](1)*_grad_test[_i][_qp](1)*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*_polar_y_grad[_qp](1)*_polar_y[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) - _grad_test[_i][_qp](1)*Utility::pow<2>(_polar_y[_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*Utility::pow<3>(_polar_y[_qp])*_test[_i][_qp]*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + (_strain[_qp](0,0)*_grad_test[_i][_qp](1)
      + _strain[_qp](2,2)*_grad_test[_i][_qp](1))*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp]
      + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (_polar_x_grad[_qp](0)*_grad_test[_i][_qp](1) + _polar_z_grad[_qp](2)*_grad_test[_i][_qp](1))*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - (_grad_test[_i][_qp](1)*Utility::pow<2>(_polar_x[_qp]) + _grad_test[_i][_qp](1)*Utility::pow<2>(_polar_z[_qp]))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*(_strain[_qp](0,0)*_polar_y[_qp]*_test[_i][_qp] + _strain[_qp](2,2)*_polar_y[_qp]*_test[_i][_qp])*(_C11[_qp]*_Q12[_qp]
      + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_polar_x_grad[_qp](0)*_polar_y[_qp]*_test[_i][_qp] + _polar_z_grad[_qp](2)*_polar_y[_qp]*_test[_i][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_test[_i][_qp]
      + _polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_test[_i][_qp])*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp] + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp]) + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp]));
  }
  else if (_component == 2)
  {
    return -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](2)*_polar_x[_qp]*_test[_i][_qp] + _polar_y_grad[_qp](2)*_polar_y[_qp]*_test[_i][_qp] + _polar_z_grad[_qp](0)*_polar_x[_qp]*_test[_i][_qp]
      + _polar_z_grad[_qp](1)*_polar_y[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](0)*_polar_x[_qp]*_polar_z[_qp] + _grad_test[_i][_qp](1)*_polar_y[_qp]*_polar_z[_qp])
      - 4*_C44[_qp]*_Q44[_qp]*(_strain[_qp](0,2)*_polar_x[_qp]*_test[_i][_qp] + _strain[_qp](1,2)*_polar_y[_qp]*_test[_i][_qp]) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(_polar_x_grad[_qp](2)*_grad_test[_i][_qp](0)
      + _polar_y_grad[_qp](2)*_grad_test[_i][_qp](1) + _polar_z_grad[_qp](0)*_grad_test[_i][_qp](0) + _polar_z_grad[_qp](1)*_grad_test[_i][_qp](1)) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(_strain[_qp](0,2)*_grad_test[_i][_qp](0)
      + _strain[_qp](1,2)*_grad_test[_i][_qp](1)) + _strain[_qp](2,2)*_grad_test[_i][_qp](2)*(_C11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_strain[_qp](2,2)*_polar_z[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]) + _polar_z_grad[_qp](2)*_grad_test[_i][_qp](2)*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*_polar_z_grad[_qp](2)*_polar_z[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) - _grad_test[_i][_qp](2)*Utility::pow<2>(_polar_z[_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*Utility::pow<3>(_polar_z[_qp])*_test[_i][_qp]*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + (_strain[_qp](0,0)*_grad_test[_i][_qp](2)
      + _strain[_qp](1,1)*_grad_test[_i][_qp](2))*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp]
      + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (_polar_x_grad[_qp](0)*_grad_test[_i][_qp](2) + _polar_y_grad[_qp](1)*_grad_test[_i][_qp](2))*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - (_grad_test[_i][_qp](2)*Utility::pow<2>(_polar_x[_qp]) + _grad_test[_i][_qp](2)*Utility::pow<2>(_polar_y[_qp]))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*(_strain[_qp](0,0)*_polar_z[_qp]*_test[_i][_qp] + _strain[_qp](1,1)*_polar_z[_qp]*_test[_i][_qp])*(_C11[_qp]*_Q12[_qp]
      + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_polar_x_grad[_qp](0)*_polar_z[_qp]*_test[_i][_qp] + _polar_y_grad[_qp](1)*_polar_z[_qp]*_test[_i][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_test[_i][_qp]
      + Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_test[_i][_qp])*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp] + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp]) + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp]));
  }
  else
    return 0.0;
}

Real
CubicParentFlexoelectricFullPDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_grad_phi[_j][_qp](1)*_polar_y[_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](2)*_polar_z[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](1)*_polar_y[_qp]*_phi[_j][_qp]
      + _grad_test[_i][_qp](2)*_polar_z[_qp]*_phi[_j][_qp]) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)
      + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)) - 2*_strain[_qp](0,0)*_phi[_j][_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 6*Utility::pow<2>(_polar_x[_qp])*_phi[_j][_qp]*_test[_i][_qp]*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) - 2*(_strain[_qp](1,1)*_phi[_j][_qp]*_test[_i][_qp]
      + _strain[_qp](2,2)*_phi[_j][_qp]*_test[_i][_qp])*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_polar_y_grad[_qp](1)*_phi[_j][_qp]*_test[_i][_qp]
      + _polar_z_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp])
      + 2*(Utility::pow<2>(_polar_y[_qp])*_phi[_j][_qp]*_test[_i][_qp] + Utility::pow<2>(_polar_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_polar_x_grad[_qp](0)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](0)*_polar_x[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](0)*_polar_x[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
  }
  else if (_component == 1)
  {
    return -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_grad_phi[_j][_qp](0)*_polar_x[_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](2)*_polar_z[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](0)*_polar_x[_qp]*_phi[_j][_qp]
      + _grad_test[_i][_qp](2)*_polar_z[_qp]*_phi[_j][_qp]) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)
      + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)) - 2*_strain[_qp](1,1)*_phi[_j][_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 6*Utility::pow<2>(_polar_y[_qp])*_phi[_j][_qp]*_test[_i][_qp]*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) - 2*(_strain[_qp](0,0)*_phi[_j][_qp]*_test[_i][_qp]
      + _strain[_qp](2,2)*_phi[_j][_qp]*_test[_i][_qp])*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_polar_x_grad[_qp](0)*_phi[_j][_qp]*_test[_i][_qp]
      + _polar_z_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp])
      + 2*(Utility::pow<2>(_polar_x[_qp])*_phi[_j][_qp]*_test[_i][_qp] + Utility::pow<2>(_polar_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_polar_y_grad[_qp](1)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](1)*_polar_y[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](1)*_polar_y[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
  }
  else if (_component == 2)
  {
    return -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_grad_phi[_j][_qp](0)*_polar_x[_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](1)*_polar_y[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](0)*_polar_x[_qp]*_phi[_j][_qp]
      + _grad_test[_i][_qp](1)*_polar_y[_qp]*_phi[_j][_qp]) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)
      + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)) - 2*_strain[_qp](2,2)*_phi[_j][_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 6*Utility::pow<2>(_polar_z[_qp])*_phi[_j][_qp]*_test[_i][_qp]*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) - 2*(_strain[_qp](0,0)*_phi[_j][_qp]*_test[_i][_qp]
      + _strain[_qp](1,1)*_phi[_j][_qp]*_test[_i][_qp])*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_polar_x_grad[_qp](0)*_phi[_j][_qp]*_test[_i][_qp]
      + _polar_y_grad[_qp](1)*_phi[_j][_qp]*_test[_i][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp])
      + 2*(Utility::pow<2>(_polar_x[_qp])*_phi[_j][_qp]*_test[_i][_qp] + Utility::pow<2>(_polar_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_polar_z_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](2)*_polar_z[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](2)*_polar_z[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
  }
  else
    return 0.0;
}

Real
CubicParentFlexoelectricFullPDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  // --- P <- u: the residual is linear in the total strain, eps -> eps + sym(e_m (x) grad phi) ---
  if (_ndisp > 0 && jvar == _disp_var[0])
    return computeQpStrainDerivative(0)*_grad_phi[_j][_qp](0) + 0.5*computeQpStrainDerivative(5)*_grad_phi[_j][_qp](1) + 0.5*computeQpStrainDerivative(4)*_grad_phi[_j][_qp](2);
  else if (_ndisp > 1 && jvar == _disp_var[1])
    return 0.5*computeQpStrainDerivative(5)*_grad_phi[_j][_qp](0) + computeQpStrainDerivative(1)*_grad_phi[_j][_qp](1) + 0.5*computeQpStrainDerivative(3)*_grad_phi[_j][_qp](2);
  else if (_ndisp > 2 && jvar == _disp_var[2])
    return 0.5*computeQpStrainDerivative(4)*_grad_phi[_j][_qp](0) + 0.5*computeQpStrainDerivative(3)*_grad_phi[_j][_qp](1) + computeQpStrainDerivative(2)*_grad_phi[_j][_qp](2);

  // --- P <- P ---
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return -4*_C44[_qp]*_strain[_qp](0,1)*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp] + _C44[_qp]*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])
      - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](1)*_phi[_j][_qp]*_test[_i][_qp] + _polar_y_grad[_qp](0)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](0)*_polar_y[_qp]*_test[_i][_qp]
      + _grad_test[_i][_qp](1)*_polar_x[_qp]*_phi[_j][_qp]) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 4*_polar_x[_qp]*_polar_y[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_grad_phi[_j][_qp](1)*_polar_x[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](0)*_polar_y[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    }
    else if (jvar == _polar_z_var)
    {
      return -4*_C44[_qp]*_strain[_qp](0,2)*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp] + _C44[_qp]*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])
      - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp] + _polar_z_grad[_qp](0)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](0)*_polar_z[_qp]*_test[_i][_qp]
      + _grad_test[_i][_qp](2)*_polar_x[_qp]*_phi[_j][_qp]) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 4*_polar_x[_qp]*_polar_z[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_grad_phi[_j][_qp](2)*_polar_x[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](0)*_polar_z[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    }
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return -4*_C44[_qp]*_strain[_qp](0,1)*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp] + _C44[_qp]*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])
      - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](1)*_phi[_j][_qp]*_test[_i][_qp] + _polar_y_grad[_qp](0)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](1)*_polar_x[_qp]*_test[_i][_qp]
      + _grad_test[_i][_qp](0)*_polar_y[_qp]*_phi[_j][_qp]) + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 4*_polar_x[_qp]*_polar_y[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_grad_phi[_j][_qp](0)*_polar_y[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](1)*_polar_x[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    }
    else if (jvar == _polar_z_var)
    {
      return -4*_C44[_qp]*_strain[_qp](1,2)*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp] + _C44[_qp]*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])
      - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_y_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp] + _polar_z_grad[_qp](1)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](1)*_polar_z[_qp]*_test[_i][_qp]
      + _grad_test[_i][_qp](2)*_polar_y[_qp]*_phi[_j][_qp]) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 4*_polar_y[_qp]*_polar_z[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_grad_phi[_j][_qp](2)*_polar_y[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](1)*_polar_z[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    }
    else
      return 0.0;
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return -4*_C44[_qp]*_strain[_qp](0,2)*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp] + _C44[_qp]*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])
      - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp] + _polar_z_grad[_qp](0)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](2)*_polar_x[_qp]*_test[_i][_qp]
      + _grad_test[_i][_qp](0)*_polar_z[_qp]*_phi[_j][_qp]) + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 4*_polar_x[_qp]*_polar_z[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_grad_phi[_j][_qp](0)*_polar_z[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](2)*_polar_x[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    }
    else if (jvar == _polar_y_var)
    {
      return -4*_C44[_qp]*_strain[_qp](1,2)*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp] + _C44[_qp]*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])
      - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_y_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp] + _polar_z_grad[_qp](1)*_phi[_j][_qp]*_test[_i][_qp] + _grad_phi[_j][_qp](2)*_polar_y[_qp]*_test[_i][_qp]
      + _grad_test[_i][_qp](1)*_polar_z[_qp]*_phi[_j][_qp]) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + 4*_polar_y[_qp]*_polar_z[_qp]*_phi[_j][_qp]*_test[_i][_qp]*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_grad_phi[_j][_qp](1)*_polar_z[_qp]*_test[_i][_qp] + _grad_test[_i][_qp](2)*_polar_y[_qp]*_phi[_j][_qp])*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
CubicParentFlexoelectricFullPDerivative::computeQpStrainDerivative(unsigned int b)
{
  if (_component == 0)
  {
    if (b == 0)
    {
      return _grad_test[_i][_qp](0)*(_C11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_x[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]);
    }
    else if (b == 1)
    {
      return _grad_test[_i][_qp](0)*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_x[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]);
    }
    else if (b == 2)
    {
      return _grad_test[_i][_qp](0)*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_x[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]);
    }
    else if (b == 3)
    {
      return 0;
    }
    else if (b == 4)
    {
      return 2*_C44[_qp]*_grad_test[_i][_qp](2)*_mu44[_qp]*_s44[_qp] - 4*_C44[_qp]*_polar_z[_qp]*_Q44[_qp]*_test[_i][_qp];
    }
    else if (b == 5)
    {
      return 2*_C44[_qp]*_grad_test[_i][_qp](1)*_mu44[_qp]*_s44[_qp] - 4*_C44[_qp]*_polar_y[_qp]*_Q44[_qp]*_test[_i][_qp];
    }
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (b == 0)
    {
      return _grad_test[_i][_qp](1)*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_y[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]);
    }
    else if (b == 1)
    {
      return _grad_test[_i][_qp](1)*(_C11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_y[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]);
    }
    else if (b == 2)
    {
      return _grad_test[_i][_qp](1)*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_y[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]);
    }
    else if (b == 3)
    {
      return 2*_C44[_qp]*_grad_test[_i][_qp](2)*_mu44[_qp]*_s44[_qp] - 4*_C44[_qp]*_polar_z[_qp]*_Q44[_qp]*_test[_i][_qp];
    }
    else if (b == 4)
    {
      return 0;
    }
    else if (b == 5)
    {
      return 2*_C44[_qp]*_grad_test[_i][_qp](0)*_mu44[_qp]*_s44[_qp] - 4*_C44[_qp]*_polar_x[_qp]*_Q44[_qp]*_test[_i][_qp];
    }
    else
      return 0.0;
  }
  else if (_component == 2)
  {
    if (b == 0)
    {
      return _grad_test[_i][_qp](2)*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_z[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]);
    }
    else if (b == 1)
    {
      return _grad_test[_i][_qp](2)*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_z[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]);
    }
    else if (b == 2)
    {
      return _grad_test[_i][_qp](2)*(_C11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - 2*_polar_z[_qp]*_test[_i][_qp]*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]);
    }
    else if (b == 3)
    {
      return 2*_C44[_qp]*_grad_test[_i][_qp](1)*_mu44[_qp]*_s44[_qp] - 4*_C44[_qp]*_polar_y[_qp]*_Q44[_qp]*_test[_i][_qp];
    }
    else if (b == 4)
    {
      return 2*_C44[_qp]*_grad_test[_i][_qp](0)*_mu44[_qp]*_s44[_qp] - 4*_C44[_qp]*_polar_x[_qp]*_Q44[_qp]*_test[_i][_qp];
    }
    else if (b == 5)
    {
      return 0;
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

void
CubicParentFlexoelectricFullPDerivative::computeOffDiagJacobianScalar(unsigned int jvar)
{
  // P <-> eps_bar. The mean strain enters the total strain like a uniform strain, so
  //     dR_i / d eps_bar_b = (d R_i / d eps) : E_b ,
  // and the eps_bar residual (integral of sigma, MOOSE's GlobalStrainUserObject) has the exact
  // transpose. MOOSE's GlobalStrain scalar kernel supplies only its own diagonal, so both
  // blocks are assembled here; each P component fills its own row, nothing twice.
  if (_global_strain_var == libMesh::invalid_uint || jvar != _global_strain_var)
    return;

  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);

  prepareMatrixTag(_assembly, _var.number(), jvar);
  for (_i = 0; _i < _test.size(); ++_i)
    for (_j = 0; _j < jv.order(); ++_j)
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpStrainDerivative(_j);
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
