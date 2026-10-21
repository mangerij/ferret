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

#include "RotatedCubicParentFlexoelectricFullPDerivative.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "Assembly.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotatedCubicParentFlexoelectricFullPDerivative);

/*
  Generated from ONE energy density (scratch_calcs/flexo_port/gen.py),

     W = 1/2 (eps - eps^Q - eps^F) : C : (eps - eps^Q - eps^F),

  evaluated on crystal-frame locals and rotated back (exact), so the residual, every Jacobian block and the energy postprocessor are exact derivatives
  of each other.  This kernel carries BOTH the electrostrictive and the flexoelectric parts
  of dW/dP: use it INSTEAD of CubicParentElasticPDerivative, never together with it.

  The flexoelectric eigenstrain is eps^F_ij = -A_ijkl dP_k/dx_l with A = s : f, where the
  flexocoupling tensor f is entered through its TENSOR components mu11 = f_1111,
  mu12 = f_1122, mu44 = f_1212 and the compliances in Voigt form (s_1212 = s44/4).
  With this sign, the published f = C : F of Gu et al., PRB 89, 174111 (2014), entered
  directly, gives their tail-to-tail Neel component.  The eigenstrain enters the
  CompositeEigenstrain with weight +1.
*/

InputParameters RotatedCubicParentFlexoelectricFullPDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual of the full elastic energy, electrostrictive AND flexoelectric, w.r.t. the polarization, for a cubic parent phase whose crystal axes are turned by the material property twist about z. Replaces CubicParentElasticPDerivative when flexoelectricity is on.");
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

RotatedCubicParentFlexoelectricFullPDerivative::RotatedCubicParentFlexoelectricFullPDerivative(const InputParameters & parameters)
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
   _twist(getMaterialProperty<Real>("twist")),
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
RotatedCubicParentFlexoelectricFullPDerivative::computeQpResidual()
{
  const Real ct = std::cos(_twist[_qp]);
  const Real st = std::sin(_twist[_qp]);
  // polarization, its gradient and the total strain in the crystal frame
  const Real p0 = _polar_x[_qp]*ct + _polar_y[_qp]*st;
  const Real p1 = -_polar_x[_qp]*st + _polar_y[_qp]*ct;
  const Real p2 = _polar_z[_qp];
  const Real g00 = _polar_x_grad[_qp](0)*Utility::pow<2>(ct) + _polar_x_grad[_qp](1)*ct*st + _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(st);
  const Real g01 = -_polar_x_grad[_qp](0)*ct*st + _polar_x_grad[_qp](1)*Utility::pow<2>(ct) - _polar_y_grad[_qp](0)*Utility::pow<2>(st) + _polar_y_grad[_qp](1)*ct*st;
  const Real g02 = _polar_x_grad[_qp](2)*ct + _polar_y_grad[_qp](2)*st;
  const Real g10 = -_polar_x_grad[_qp](0)*ct*st - _polar_x_grad[_qp](1)*Utility::pow<2>(st) + _polar_y_grad[_qp](0)*Utility::pow<2>(ct) + _polar_y_grad[_qp](1)*ct*st;
  const Real g11 = _polar_x_grad[_qp](0)*Utility::pow<2>(st) - _polar_x_grad[_qp](1)*ct*st - _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(ct);
  const Real g12 = -_polar_x_grad[_qp](2)*st + _polar_y_grad[_qp](2)*ct;
  const Real g20 = _polar_z_grad[_qp](0)*ct + _polar_z_grad[_qp](1)*st;
  const Real g21 = -_polar_z_grad[_qp](0)*st + _polar_z_grad[_qp](1)*ct;
  const Real g22 = _polar_z_grad[_qp](2);
  const Real e00 = _strain[_qp](0,0)*Utility::pow<2>(ct) + 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(st);
  const Real e01 = -_strain[_qp](0,0)*ct*st + _strain[_qp](0,1)*Utility::pow<2>(ct) - _strain[_qp](0,1)*Utility::pow<2>(st) + _strain[_qp](1,1)*ct*st;
  const Real e02 = _strain[_qp](0,2)*ct + _strain[_qp](1,2)*st;
  const Real e11 = _strain[_qp](0,0)*Utility::pow<2>(st) - 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(ct);
  const Real e12 = -_strain[_qp](0,2)*st + _strain[_qp](1,2)*ct;
  const Real e22 = _strain[_qp](2,2);
  // test-function gradient in the crystal frame
  const Real gt0 = _grad_test[_i][_qp](0)*ct + _grad_test[_i][_qp](1)*st;
  const Real gt1 = -_grad_test[_i][_qp](0)*st + _grad_test[_i][_qp](1)*ct;
  const Real gt2 = _grad_test[_i][_qp](2);

  if (_component == 2)
    return -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_test[_i][_qp]*g02*p0 + _test[_i][_qp]*g12*p1 + _test[_i][_qp]*g20*p0 + _test[_i][_qp]*g21*p1 + gt0*p0*p2 + gt1*p1*p2)
      - 4*_C44[_qp]*_Q44[_qp]*(_test[_i][_qp]*e02*p0 + _test[_i][_qp]*e12*p1) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(g02*gt0 + g12*gt1
      + g20*gt0 + g21*gt1) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(e02*gt0 + e12*gt1) - 2*_test[_i][_qp]*e22*p2*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) - 2*_test[_i][_qp]*g22*p2*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*_test[_i][_qp]*Utility::pow<3>(p2)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + e22*gt2*(_C11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + g22*gt2*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - gt2*Utility::pow<2>(p2)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + (e00*gt2 + e11*gt2)*(_C11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp]
      + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (g00*gt2 + g11*gt2)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - (gt2*Utility::pow<2>(p0) + gt2*Utility::pow<2>(p1))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp])
      - 2*(_test[_i][_qp]*e00*p2 + _test[_i][_qp]*e11*p2)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_test[_i][_qp]*g00*p2 + _test[_i][_qp]*g11*p2)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_test[_i][_qp]*Utility::pow<2>(p0)*p2 + _test[_i][_qp]*Utility::pow<2>(p1)*p2)*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp]));

  // x and y mix under the twist: (sim) = R(twist) (crystal)
  const Real r0 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_test[_i][_qp]*g01*p1 + _test[_i][_qp]*g02*p2 + _test[_i][_qp]*g10*p1 + _test[_i][_qp]*g20*p2 + gt1*p0*p1 + gt2*p0*p2)
      - 4*_C44[_qp]*_Q44[_qp]*(_test[_i][_qp]*e01*p1 + _test[_i][_qp]*e02*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(g01*gt1 + g02*gt2
      + g10*gt1 + g20*gt2) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(e01*gt1 + e02*gt2) - 2*_test[_i][_qp]*e00*p0*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) - 2*_test[_i][_qp]*g00*p0*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*_test[_i][_qp]*Utility::pow<3>(p0)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + e00*gt0*(_C11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + g00*gt0*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - gt0*Utility::pow<2>(p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + (e11*gt0 + e22*gt0)*(_C11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp]
      + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (g11*gt0 + g22*gt0)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - (gt0*Utility::pow<2>(p1) + gt0*Utility::pow<2>(p2))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp])
      - 2*(_test[_i][_qp]*e11*p0 + _test[_i][_qp]*e22*p0)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_test[_i][_qp]*g11*p0 + _test[_i][_qp]*g22*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_test[_i][_qp]*p0*Utility::pow<2>(p1) + _test[_i][_qp]*p0*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp]));
  const Real r1 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_test[_i][_qp]*g01*p0 + _test[_i][_qp]*g10*p0 + _test[_i][_qp]*g12*p2 + _test[_i][_qp]*g21*p2 + gt0*p0*p1 + gt2*p1*p2)
      - 4*_C44[_qp]*_Q44[_qp]*(_test[_i][_qp]*e01*p0 + _test[_i][_qp]*e12*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(g01*gt0 + g10*gt0
      + g12*gt2 + g21*gt2) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(e01*gt0 + e12*gt2) - 2*_test[_i][_qp]*e11*p1*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) - 2*_test[_i][_qp]*g11*p1*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*_test[_i][_qp]*Utility::pow<3>(p1)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + e11*gt1*(_C11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + g11*gt1*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - gt1*Utility::pow<2>(p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + (e00*gt1 + e22*gt1)*(_C11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp]
      + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (g00*gt1 + g22*gt1)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - (gt1*Utility::pow<2>(p0) + gt1*Utility::pow<2>(p2))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp])
      - 2*(_test[_i][_qp]*e00*p1 + _test[_i][_qp]*e22*p1)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_test[_i][_qp]*g00*p1 + _test[_i][_qp]*g22*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_test[_i][_qp]*Utility::pow<2>(p0)*p1 + _test[_i][_qp]*p1*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp]));
  if (_component == 0)
    return ct*r0 - st*r1;
  else
    return st*r0 + ct*r1;
}

Real
RotatedCubicParentFlexoelectricFullPDerivative::computeQpJacobian()
{
  const Real ct = std::cos(_twist[_qp]);
  const Real st = std::sin(_twist[_qp]);
  // polarization, its gradient and the total strain in the crystal frame
  const Real p0 = _polar_x[_qp]*ct + _polar_y[_qp]*st;
  const Real p1 = -_polar_x[_qp]*st + _polar_y[_qp]*ct;
  const Real p2 = _polar_z[_qp];
  const Real g00 = _polar_x_grad[_qp](0)*Utility::pow<2>(ct) + _polar_x_grad[_qp](1)*ct*st + _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(st);
  const Real g01 = -_polar_x_grad[_qp](0)*ct*st + _polar_x_grad[_qp](1)*Utility::pow<2>(ct) - _polar_y_grad[_qp](0)*Utility::pow<2>(st) + _polar_y_grad[_qp](1)*ct*st;
  const Real g02 = _polar_x_grad[_qp](2)*ct + _polar_y_grad[_qp](2)*st;
  const Real g10 = -_polar_x_grad[_qp](0)*ct*st - _polar_x_grad[_qp](1)*Utility::pow<2>(st) + _polar_y_grad[_qp](0)*Utility::pow<2>(ct) + _polar_y_grad[_qp](1)*ct*st;
  const Real g11 = _polar_x_grad[_qp](0)*Utility::pow<2>(st) - _polar_x_grad[_qp](1)*ct*st - _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(ct);
  const Real g12 = -_polar_x_grad[_qp](2)*st + _polar_y_grad[_qp](2)*ct;
  const Real g20 = _polar_z_grad[_qp](0)*ct + _polar_z_grad[_qp](1)*st;
  const Real g21 = -_polar_z_grad[_qp](0)*st + _polar_z_grad[_qp](1)*ct;
  const Real g22 = _polar_z_grad[_qp](2);
  const Real e00 = _strain[_qp](0,0)*Utility::pow<2>(ct) + 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(st);
  const Real e01 = -_strain[_qp](0,0)*ct*st + _strain[_qp](0,1)*Utility::pow<2>(ct) - _strain[_qp](0,1)*Utility::pow<2>(st) + _strain[_qp](1,1)*ct*st;
  const Real e02 = _strain[_qp](0,2)*ct + _strain[_qp](1,2)*st;
  const Real e11 = _strain[_qp](0,0)*Utility::pow<2>(st) - 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(ct);
  const Real e12 = -_strain[_qp](0,2)*st + _strain[_qp](1,2)*ct;
  const Real e22 = _strain[_qp](2,2);
  // test-function gradient in the crystal frame
  const Real gt0 = _grad_test[_i][_qp](0)*ct + _grad_test[_i][_qp](1)*st;
  const Real gt1 = -_grad_test[_i][_qp](0)*st + _grad_test[_i][_qp](1)*ct;
  const Real gt2 = _grad_test[_i][_qp](2);
  // shape-function gradient in the crystal frame
  const Real gp0 = _grad_phi[_j][_qp](0)*ct + _grad_phi[_j][_qp](1)*st;
  const Real gp1 = -_grad_phi[_j][_qp](0)*st + _grad_phi[_j][_qp](1)*ct;
  const Real gp2 = _grad_phi[_j][_qp](2);

  if (_component == 0)
  {
    const Real j00 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt1*p1 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp1*p1 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp1*gt1
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e00*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp0*gt0*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e11
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*gt0*p0 + _test[_i][_qp]*gp0*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    const Real j01 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1) + _C44[_qp]*gp0*gt1*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp1*gt0*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    const Real j10 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0) + _C44[_qp]*gp1*gt0*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp0*gt1*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    const Real j11 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt0*p0 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp0*p0 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp0*gt0
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e11*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp1*gt1*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e00
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*gt1*p1 + _test[_i][_qp]*gp1*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    return Utility::pow<2>(ct)*j00 - ct*j01*st - ct*j10*st + j11*Utility::pow<2>(st);
  }
  else if (_component == 1)
  {
    const Real j00 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt1*p1 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp1*p1 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp1*gt1
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e00*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp0*gt0*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e11
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*gt0*p0 + _test[_i][_qp]*gp0*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    const Real j01 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1) + _C44[_qp]*gp0*gt1*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp1*gt0*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    const Real j10 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0) + _C44[_qp]*gp1*gt0*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp0*gt1*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    const Real j11 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt0*p0 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp0*p0 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp0*gt0
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e11*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp1*gt1*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e00
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*gt1*p1 + _test[_i][_qp]*gp1*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    return Utility::pow<2>(ct)*j11 + ct*j01*st + ct*j10*st + j00*Utility::pow<2>(st);
  }
  else if (_component == 2)
  {
    const Real j22 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt0*p0 + _phi[_j][_qp]*gt1*p1 + _test[_i][_qp]*gp0*p0 + _test[_i][_qp]*gp1*p1) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp0*gt0
      + gp1*gt1) - 2*_phi[_j][_qp]*_test[_i][_qp]*e22*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp2*gt2*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e00
      + _phi[_j][_qp]*_test[_i][_qp]*e11)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*_test[_i][_qp]*g11)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g22 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp2*p2)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
    return j22;
  }
  else
    return 0.0;
}

Real
RotatedCubicParentFlexoelectricFullPDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  // --- P <- u: the residual is linear in the total strain, eps -> eps + sym(e_m (x) grad phi) ---
  if (_ndisp > 0 && jvar == _disp_var[0])
    return computeQpStrainDerivative(0)*_grad_phi[_j][_qp](0) + 0.5*computeQpStrainDerivative(5)*_grad_phi[_j][_qp](1) + 0.5*computeQpStrainDerivative(4)*_grad_phi[_j][_qp](2);
  else if (_ndisp > 1 && jvar == _disp_var[1])
    return 0.5*computeQpStrainDerivative(5)*_grad_phi[_j][_qp](0) + computeQpStrainDerivative(1)*_grad_phi[_j][_qp](1) + 0.5*computeQpStrainDerivative(3)*_grad_phi[_j][_qp](2);
  else if (_ndisp > 2 && jvar == _disp_var[2])
    return 0.5*computeQpStrainDerivative(4)*_grad_phi[_j][_qp](0) + 0.5*computeQpStrainDerivative(3)*_grad_phi[_j][_qp](1) + computeQpStrainDerivative(2)*_grad_phi[_j][_qp](2);

  // --- P <- P ---
  const Real ct = std::cos(_twist[_qp]);
  const Real st = std::sin(_twist[_qp]);
  // polarization, its gradient and the total strain in the crystal frame
  const Real p0 = _polar_x[_qp]*ct + _polar_y[_qp]*st;
  const Real p1 = -_polar_x[_qp]*st + _polar_y[_qp]*ct;
  const Real p2 = _polar_z[_qp];
  const Real g00 = _polar_x_grad[_qp](0)*Utility::pow<2>(ct) + _polar_x_grad[_qp](1)*ct*st + _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(st);
  const Real g01 = -_polar_x_grad[_qp](0)*ct*st + _polar_x_grad[_qp](1)*Utility::pow<2>(ct) - _polar_y_grad[_qp](0)*Utility::pow<2>(st) + _polar_y_grad[_qp](1)*ct*st;
  const Real g02 = _polar_x_grad[_qp](2)*ct + _polar_y_grad[_qp](2)*st;
  const Real g10 = -_polar_x_grad[_qp](0)*ct*st - _polar_x_grad[_qp](1)*Utility::pow<2>(st) + _polar_y_grad[_qp](0)*Utility::pow<2>(ct) + _polar_y_grad[_qp](1)*ct*st;
  const Real g11 = _polar_x_grad[_qp](0)*Utility::pow<2>(st) - _polar_x_grad[_qp](1)*ct*st - _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(ct);
  const Real g12 = -_polar_x_grad[_qp](2)*st + _polar_y_grad[_qp](2)*ct;
  const Real g20 = _polar_z_grad[_qp](0)*ct + _polar_z_grad[_qp](1)*st;
  const Real g21 = -_polar_z_grad[_qp](0)*st + _polar_z_grad[_qp](1)*ct;
  const Real g22 = _polar_z_grad[_qp](2);
  const Real e00 = _strain[_qp](0,0)*Utility::pow<2>(ct) + 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(st);
  const Real e01 = -_strain[_qp](0,0)*ct*st + _strain[_qp](0,1)*Utility::pow<2>(ct) - _strain[_qp](0,1)*Utility::pow<2>(st) + _strain[_qp](1,1)*ct*st;
  const Real e02 = _strain[_qp](0,2)*ct + _strain[_qp](1,2)*st;
  const Real e11 = _strain[_qp](0,0)*Utility::pow<2>(st) - 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(ct);
  const Real e12 = -_strain[_qp](0,2)*st + _strain[_qp](1,2)*ct;
  const Real e22 = _strain[_qp](2,2);
  // test-function gradient in the crystal frame
  const Real gt0 = _grad_test[_i][_qp](0)*ct + _grad_test[_i][_qp](1)*st;
  const Real gt1 = -_grad_test[_i][_qp](0)*st + _grad_test[_i][_qp](1)*ct;
  const Real gt2 = _grad_test[_i][_qp](2);
  // shape-function gradient in the crystal frame
  const Real gp0 = _grad_phi[_j][_qp](0)*ct + _grad_phi[_j][_qp](1)*st;
  const Real gp1 = -_grad_phi[_j][_qp](0)*st + _grad_phi[_j][_qp](1)*ct;
  const Real gp2 = _grad_phi[_j][_qp](2);

  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      const Real j00 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt1*p1 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp1*p1 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp1*gt1
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e00*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp0*gt0*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e11
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*gt0*p0 + _test[_i][_qp]*gp0*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j01 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1) + _C44[_qp]*gp0*gt1*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp1*gt0*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j10 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0) + _C44[_qp]*gp1*gt0*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp0*gt1*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j11 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt0*p0 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp0*p0 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp0*gt0
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e11*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp1*gt1*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e00
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*gt1*p1 + _test[_i][_qp]*gp1*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      return Utility::pow<2>(ct)*j01 + ct*j00*st - ct*j11*st - j10*Utility::pow<2>(st);
    }
    else if (jvar == _polar_z_var)
    {
      const Real j02 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e02 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g02 + _phi[_j][_qp]*_test[_i][_qp]*g20
      + _phi[_j][_qp]*gt2*p0 + _test[_i][_qp]*gp0*p2) + _C44[_qp]*gp0*gt2*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp2*gt0*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt0*p2 + _test[_i][_qp]*gp2*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j12 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e12 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g12 + _phi[_j][_qp]*_test[_i][_qp]*g21
      + _phi[_j][_qp]*gt2*p1 + _test[_i][_qp]*gp1*p2) + _C44[_qp]*gp1*gt2*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p1*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp2*gt1*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt1*p2 + _test[_i][_qp]*gp2*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      return ct*j02 - j12*st;
    }
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      const Real j00 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt1*p1 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp1*p1 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp1*gt1
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e00*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp0*gt0*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e11
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*gt0*p0 + _test[_i][_qp]*gp0*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j01 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1) + _C44[_qp]*gp0*gt1*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp1*gt0*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j10 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e01 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g01 + _phi[_j][_qp]*_test[_i][_qp]*g10
      + _phi[_j][_qp]*gt0*p1 + _test[_i][_qp]*gp1*p0) + _C44[_qp]*gp1*gt0*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p1*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp0*gt1*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt1*p0 + _test[_i][_qp]*gp0*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j11 = -2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*gt0*p0 + _phi[_j][_qp]*gt2*p2 + _test[_i][_qp]*gp0*p0 + _test[_i][_qp]*gp2*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(gp0*gt0
      + gp2*gt2) - 2*_phi[_j][_qp]*_test[_i][_qp]*e11*(_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp]) + 6*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p1)*(_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 4*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_Q12[_qp])) + gp1*gt1*(_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 8*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 8*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*e00
      + _phi[_j][_qp]*_test[_i][_qp]*e22)*(_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp]) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g00 + _phi[_j][_qp]*_test[_i][_qp]*g22)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + 2*(_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p0) + _phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) - 2*(_phi[_j][_qp]*_test[_i][_qp]*g11 + _phi[_j][_qp]*gt1*p1 + _test[_i][_qp]*gp1*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      return Utility::pow<2>(ct)*j10 + ct*j00*st - ct*j11*st - j01*Utility::pow<2>(st);
    }
    else if (jvar == _polar_z_var)
    {
      const Real j02 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e02 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g02 + _phi[_j][_qp]*_test[_i][_qp]*g20
      + _phi[_j][_qp]*gt2*p0 + _test[_i][_qp]*gp0*p2) + _C44[_qp]*gp0*gt2*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp2*gt0*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt0*p2 + _test[_i][_qp]*gp2*p0)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j12 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e12 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g12 + _phi[_j][_qp]*_test[_i][_qp]*g21
      + _phi[_j][_qp]*gt2*p1 + _test[_i][_qp]*gp1*p2) + _C44[_qp]*gp1*gt2*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p1*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp2*gt1*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt1*p2 + _test[_i][_qp]*gp2*p1)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      return ct*j12 + j02*st;
    }
    else
      return 0.0;
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      const Real j20 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e02 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g02 + _phi[_j][_qp]*_test[_i][_qp]*g20
      + _phi[_j][_qp]*gt0*p2 + _test[_i][_qp]*gp2*p0) + _C44[_qp]*gp2*gt0*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp0*gt2*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt2*p0 + _test[_i][_qp]*gp0*p2)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j21 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e12 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g12 + _phi[_j][_qp]*_test[_i][_qp]*g21
      + _phi[_j][_qp]*gt1*p2 + _test[_i][_qp]*gp2*p1) + _C44[_qp]*gp2*gt1*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p1*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp1*gt2*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt2*p1 + _test[_i][_qp]*gp1*p2)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      return ct*j20 - j21*st;
    }
    else if (jvar == _polar_y_var)
    {
      const Real j20 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e02 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g02 + _phi[_j][_qp]*_test[_i][_qp]*g20
      + _phi[_j][_qp]*gt0*p2 + _test[_i][_qp]*gp2*p0) + _C44[_qp]*gp2*gt0*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p0*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp0*gt2*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt2*p0 + _test[_i][_qp]*gp0*p2)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      const Real j21 = -4*_C44[_qp]*_phi[_j][_qp]*_Q44[_qp]*_test[_i][_qp]*e12 - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_phi[_j][_qp]*_test[_i][_qp]*g12 + _phi[_j][_qp]*_test[_i][_qp]*g21
      + _phi[_j][_qp]*gt1*p2 + _test[_i][_qp]*gp2*p1) + _C44[_qp]*gp2*gt1*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp]) + 4*_phi[_j][_qp]*_test[_i][_qp]*p1*p2*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + gp1*gt2*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]))
      - 2*(_phi[_j][_qp]*gt2*p1 + _test[_i][_qp]*gp1*p2)*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp]
      + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]);
      return ct*j21 + j20*st;
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
RotatedCubicParentFlexoelectricFullPDerivative::computeQpStrainDerivative(unsigned int b)
{
  const Real ct = std::cos(_twist[_qp]);
  const Real st = std::sin(_twist[_qp]);
  // polarization, its gradient and the total strain in the crystal frame
  const Real p0 = _polar_x[_qp]*ct + _polar_y[_qp]*st;
  const Real p1 = -_polar_x[_qp]*st + _polar_y[_qp]*ct;
  const Real p2 = _polar_z[_qp];
  const Real g00 = _polar_x_grad[_qp](0)*Utility::pow<2>(ct) + _polar_x_grad[_qp](1)*ct*st + _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(st);
  const Real g01 = -_polar_x_grad[_qp](0)*ct*st + _polar_x_grad[_qp](1)*Utility::pow<2>(ct) - _polar_y_grad[_qp](0)*Utility::pow<2>(st) + _polar_y_grad[_qp](1)*ct*st;
  const Real g02 = _polar_x_grad[_qp](2)*ct + _polar_y_grad[_qp](2)*st;
  const Real g10 = -_polar_x_grad[_qp](0)*ct*st - _polar_x_grad[_qp](1)*Utility::pow<2>(st) + _polar_y_grad[_qp](0)*Utility::pow<2>(ct) + _polar_y_grad[_qp](1)*ct*st;
  const Real g11 = _polar_x_grad[_qp](0)*Utility::pow<2>(st) - _polar_x_grad[_qp](1)*ct*st - _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(ct);
  const Real g12 = -_polar_x_grad[_qp](2)*st + _polar_y_grad[_qp](2)*ct;
  const Real g20 = _polar_z_grad[_qp](0)*ct + _polar_z_grad[_qp](1)*st;
  const Real g21 = -_polar_z_grad[_qp](0)*st + _polar_z_grad[_qp](1)*ct;
  const Real g22 = _polar_z_grad[_qp](2);
  // test-function gradient in the crystal frame
  const Real gt0 = _grad_test[_i][_qp](0)*ct + _grad_test[_i][_qp](1)*st;
  const Real gt1 = -_grad_test[_i][_qp](0)*st + _grad_test[_i][_qp](1)*ct;
  const Real gt2 = _grad_test[_i][_qp](2);
  // the unit mean strain E_b (MOOSE ordering xx yy zz yz xz xy), turned into the crystal frame
  const Real b00 = (b == 0 ? 1.0 : 0.0);
  const Real b11 = (b == 1 ? 1.0 : 0.0);
  const Real b22 = (b == 2 ? 1.0 : 0.0);
  const Real b12 = (b == 3 ? 1.0 : 0.0);
  const Real b02 = (b == 4 ? 1.0 : 0.0);
  const Real b01 = (b == 5 ? 1.0 : 0.0);
  const Real u00 = b00*Utility::pow<2>(ct) + 2*b01*ct*st + b11*Utility::pow<2>(st);
  const Real u01 = -b00*ct*st + b01*Utility::pow<2>(ct) - b01*Utility::pow<2>(st) + b11*ct*st;
  const Real u02 = b02*ct + b12*st;
  const Real u11 = b00*Utility::pow<2>(st) - 2*b01*ct*st + b11*Utility::pow<2>(ct);
  const Real u12 = -b02*st + b12*ct;
  const Real u22 = b22;

  if (_component == 2)
    return -4*_C44[_qp]*_Q44[_qp]*_test[_i][_qp]*p0*u02 - 4*_C44[_qp]*_Q44[_qp]*_test[_i][_qp]*p1*u12 + 2*_C44[_qp]*gt0*_mu44[_qp]*_s44[_qp]*u02 + 2*_C44[_qp]*gt1*_mu44[_qp]*_s44[_qp]*u12
      - 2*_test[_i][_qp]*p2*(_C11[_qp]*_Q11[_qp]*u22 + _C11[_qp]*_Q12[_qp]*u00 + _C11[_qp]*_Q12[_qp]*u11 + _C12[_qp]*_Q11[_qp]*u00 + _C12[_qp]*_Q11[_qp]*u11 + _C12[_qp]*_Q12[_qp]*u00
      + _C12[_qp]*_Q12[_qp]*u11 + 2*_C12[_qp]*_Q12[_qp]*u22) + gt2*(_C11[_qp]*_mu11[_qp]*_s11[_qp]*u22 + _C11[_qp]*_mu11[_qp]*_s12[_qp]*u00 + _C11[_qp]*_mu11[_qp]*_s12[_qp]*u11
      + _C11[_qp]*_mu12[_qp]*_s11[_qp]*u00 + _C11[_qp]*_mu12[_qp]*_s11[_qp]*u11 + _C11[_qp]*_mu12[_qp]*_s12[_qp]*u00 + _C11[_qp]*_mu12[_qp]*_s12[_qp]*u11 + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp]*u22
      + _C12[_qp]*_mu11[_qp]*_s11[_qp]*u00 + _C12[_qp]*_mu11[_qp]*_s11[_qp]*u11 + _C12[_qp]*_mu11[_qp]*_s12[_qp]*u00 + _C12[_qp]*_mu11[_qp]*_s12[_qp]*u11 + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp]*u22
      + _C12[_qp]*_mu12[_qp]*_s11[_qp]*u00 + _C12[_qp]*_mu12[_qp]*_s11[_qp]*u11 + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp]*u22 + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u00 + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u11
      + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u22);

  // x and y mix under the twist: (sim) = R(twist) (crystal)
  const Real r0 = -4*_C44[_qp]*_Q44[_qp]*_test[_i][_qp]*p1*u01 - 4*_C44[_qp]*_Q44[_qp]*_test[_i][_qp]*p2*u02 + 2*_C44[_qp]*gt1*_mu44[_qp]*_s44[_qp]*u01 + 2*_C44[_qp]*gt2*_mu44[_qp]*_s44[_qp]*u02
      - 2*_test[_i][_qp]*p0*(_C11[_qp]*_Q11[_qp]*u00 + _C11[_qp]*_Q12[_qp]*u11 + _C11[_qp]*_Q12[_qp]*u22 + _C12[_qp]*_Q11[_qp]*u11 + _C12[_qp]*_Q11[_qp]*u22 + 2*_C12[_qp]*_Q12[_qp]*u00
      + _C12[_qp]*_Q12[_qp]*u11 + _C12[_qp]*_Q12[_qp]*u22) + gt0*(_C11[_qp]*_mu11[_qp]*_s11[_qp]*u00 + _C11[_qp]*_mu11[_qp]*_s12[_qp]*u11 + _C11[_qp]*_mu11[_qp]*_s12[_qp]*u22
      + _C11[_qp]*_mu12[_qp]*_s11[_qp]*u11 + _C11[_qp]*_mu12[_qp]*_s11[_qp]*u22 + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp]*u00 + _C11[_qp]*_mu12[_qp]*_s12[_qp]*u11 + _C11[_qp]*_mu12[_qp]*_s12[_qp]*u22
      + _C12[_qp]*_mu11[_qp]*_s11[_qp]*u11 + _C12[_qp]*_mu11[_qp]*_s11[_qp]*u22 + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp]*u00 + _C12[_qp]*_mu11[_qp]*_s12[_qp]*u11 + _C12[_qp]*_mu11[_qp]*_s12[_qp]*u22
      + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp]*u00 + _C12[_qp]*_mu12[_qp]*_s11[_qp]*u11 + _C12[_qp]*_mu12[_qp]*_s11[_qp]*u22 + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u00 + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u11
      + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u22);
  const Real r1 = -4*_C44[_qp]*_Q44[_qp]*_test[_i][_qp]*p0*u01 - 4*_C44[_qp]*_Q44[_qp]*_test[_i][_qp]*p2*u12 + 2*_C44[_qp]*gt0*_mu44[_qp]*_s44[_qp]*u01 + 2*_C44[_qp]*gt2*_mu44[_qp]*_s44[_qp]*u12
      - 2*_test[_i][_qp]*p1*(_C11[_qp]*_Q11[_qp]*u11 + _C11[_qp]*_Q12[_qp]*u00 + _C11[_qp]*_Q12[_qp]*u22 + _C12[_qp]*_Q11[_qp]*u00 + _C12[_qp]*_Q11[_qp]*u22 + _C12[_qp]*_Q12[_qp]*u00
      + 2*_C12[_qp]*_Q12[_qp]*u11 + _C12[_qp]*_Q12[_qp]*u22) + gt1*(_C11[_qp]*_mu11[_qp]*_s11[_qp]*u11 + _C11[_qp]*_mu11[_qp]*_s12[_qp]*u00 + _C11[_qp]*_mu11[_qp]*_s12[_qp]*u22
      + _C11[_qp]*_mu12[_qp]*_s11[_qp]*u00 + _C11[_qp]*_mu12[_qp]*_s11[_qp]*u22 + _C11[_qp]*_mu12[_qp]*_s12[_qp]*u00 + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp]*u11 + _C11[_qp]*_mu12[_qp]*_s12[_qp]*u22
      + _C12[_qp]*_mu11[_qp]*_s11[_qp]*u00 + _C12[_qp]*_mu11[_qp]*_s11[_qp]*u22 + _C12[_qp]*_mu11[_qp]*_s12[_qp]*u00 + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp]*u11 + _C12[_qp]*_mu11[_qp]*_s12[_qp]*u22
      + _C12[_qp]*_mu12[_qp]*_s11[_qp]*u00 + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp]*u11 + _C12[_qp]*_mu12[_qp]*_s11[_qp]*u22 + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u00 + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u11
      + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]*u22);
  if (_component == 0)
    return ct*r0 - st*r1;
  else
    return st*r0 + ct*r1;
}

void
RotatedCubicParentFlexoelectricFullPDerivative::computeOffDiagJacobianScalar(unsigned int jvar)
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
