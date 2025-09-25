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

#include "PolarElectricPStrongEConst.h"

class PolarElectricPStrongEConst;

registerMooseObject("FerretApp", PolarElectricPStrongEConst);

InputParameters PolarElectricPStrongEConst::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to -P$*$E term in the total energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("E_x", "The fixed x-component of the electric field variable");
  params.addRequiredCoupledVar("E_y", "The fixed y-component of the electric field variable");
  params.addRequiredCoupledVar("E_z", "The fixed z-component of the electric field variable");
  return params;
}

PolarElectricPStrongEConst::PolarElectricPStrongEConst(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _E_x(coupledValue("E_x")),
   _E_y(coupledValue("E_y")),
   _E_z(coupledValue("E_z"))
{
}

Real
PolarElectricPStrongEConst::computeQpResidual()
{
    Real RpolarP = 0.0;
    RealVectorValue Evec(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
    RpolarP += -(Evec(_component)) * _test[_i][_qp];

    ///  Moose::out << "\n R_polarP-"; std::cout << _component << " = " << RpolarP;

    return RpolarP;
}

Real
PolarElectricPStrongEConst::computeQpJacobian()
{
  return 0.0;
}
