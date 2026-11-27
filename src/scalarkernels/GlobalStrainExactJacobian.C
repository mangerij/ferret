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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "GlobalStrainExactJacobian.h"
#include "GlobalStrainUserObjectInterface.h"

registerMooseObject("FerretApp", GlobalStrainExactJacobian);

InputParameters
GlobalStrainExactJacobian::validParams()
{
  InputParameters params = GlobalStrain::validParams();
  params.addClassDescription("MOOSE's GlobalStrain scalar kernel with the exact eps_bar <- eps_bar "
                             "Jacobian (MOOSE's shear columns are half the true derivative). Use "
                             "it in place of the GlobalStrain action's scalar kernel.");
  return params;
}

GlobalStrainExactJacobian::GlobalStrainExactJacobian(const InputParameters & parameters)
  : GlobalStrain(parameters)
{
}

void
GlobalStrainExactJacobian::computeJacobian()
{
  prepareMatrixTag(_assembly, _var.number(), _var.number());
  for (_i = 0; _i < _local_ke.m(); ++_i)
    for (_j = 0; _j < _local_ke.m(); ++_j)
    {
      const unsigned int a = _components[_i].first, b = _components[_i].second;
      const unsigned int k = _components[_j].first, l = _components[_j].second;
      // as MOOSE, no periodic-direction check here (avoids a zero pivot)
      _local_ke(_i, _j) += _pst_jacobian(a, b, k, l);
      if (k != l)
        _local_ke(_i, _j) += _pst_jacobian(a, b, l, k);
    }
  accumulateTaggedLocalMatrix();
}
