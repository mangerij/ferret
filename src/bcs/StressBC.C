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

#include "StressBC.h"

registerMooseObject("FerretApp", StressBC);

InputParameters StressBC::validParams()
{
    InputParameters params = IntegratedBC::validParams();
    params.addRequiredParam<int>("component","Which component(0 for x, 1 for y, 2 for z) in traction is used");
    params.addParam<std::vector<Real> >("boundary_stress", "Boundary stress: s11, s22, s33, s23, s13, s12");
    params.addCoupledVar("boundary_stress_vars", "Variable names for the: s11, s22, s33");
    params.addParam<bool>("convert_to_gpa", false, "Convert the input amounts to GPa from Pa");
    params.addParam<Real>("prefactStress", 1.0, "A prefactor multiplier");
    return params;
}

StressBC::StressBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _stress_vector(getParam<std::vector<Real> >("boundary_stress")),
    _Jacobian_mult(getMaterialProperty<RankFourTensor>("Jacobian_mult")),
    _component(getParam<int>("component")),
    _convert_to_gpa(getParam<bool>("convert_to_gpa")),
    _prefactStress(getParam<Real>("prefactStress")),
    _multiplier(1)
{
  if(isParamValid("boundary_stress"))
    _boundary_stress.fillFromInputVector(_stress_vector);
  else if(isCoupled("boundary_stress_vars"))
  {
    int n = coupledComponents("boundary_stress_vars");

    if(n > 3)
      mooseError("Can only take the diagonals as coupled values right now!");

    _boundary_stress_vars.resize(n);

    for (unsigned int i=0; i<_boundary_stress_vars.size(); ++i)
      _boundary_stress_vars[i] = &coupledValue("boundary_stress_vars", i);
  }
  else
    mooseError("Must provide something to StressBC!");

  if(_convert_to_gpa)
    _multiplier = 1e-9;
}

Real
StressBC::computeQpResidual()
{
  /// If nothing was coupled this will be a no-op
  for(unsigned int i=0; i<_boundary_stress_vars.size(); i++)
    _boundary_stress(i, i) = _multiplier * (*_boundary_stress_vars[i])[_qp];

  return -_prefactStress * _test[_i][_qp] * (_boundary_stress.row(_component)*_normals[_qp]);
}

Real
StressBC::computeQpJacobian()
{
  return 0;
}
