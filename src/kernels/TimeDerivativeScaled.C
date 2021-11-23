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

#include "TimeDerivativeScaled.h"

#include "MooseVariable.h"

registerMooseObject("FerretApp", TimeDerivativeScaled);

InputParameters TimeDerivativeScaled::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  params.addParam<Real>("time_scale",1.0,"the time_scale of the unit");
  return params;
}

TimeDerivativeScaled::TimeDerivativeScaled(const InputParameters & parameters) :
    TimeKernel(parameters),
    _lumping(getParam<bool>("lumping")),
    _time_scale(getParam<Real>("time_scale"))
{
}

Real
TimeDerivativeScaled::computeQpResidual()
{
  return _time_scale * _test[_i][_qp] * _u_dot[_qp];
}

Real
TimeDerivativeScaled::computeQpJacobian()
{
  return _time_scale * _test[_i][_qp]*_phi[_j][_qp] * _du_dot_du[_qp];
}

void
TimeDerivativeScaled::computeJacobian()
{
  if (_lumping)
  {
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
          ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
        }
  }
  else
    TimeKernel::computeJacobian();
}
