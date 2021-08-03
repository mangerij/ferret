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

   You should have received a co_antiferrodis_A_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "CorrectionFactorDistortBFO.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CorrectionFactorDistortBFO);

template<>
InputParameters validParams<CorrectionFactorDistortBFO>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a correction factor for the local free energy");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodis_Aization");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodis_Aization");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodis_Aization");
  return params;
}

CorrectionFactorDistortBFO::CorrectionFactorDistortBFO(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _f1(getMaterialProperty<Real>("f1"))
{
}

Real
CorrectionFactorDistortBFO::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (_f1[_qp]*(-4*Utility::pow<3>(_antiferrodis_A_x[_qp]) + 8*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 8*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (_f1[_qp]*(8*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp] - 4*Utility::pow<3>(_antiferrodis_A_y[_qp]) + 8*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (_f1[_qp]*(8*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp] + 8*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp] - 4*Utility::pow<3>(_antiferrodis_A_z[_qp])));
  }
  else
    return 0.0;
}

Real
CorrectionFactorDistortBFO::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (_f1[_qp]*(-12*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 8*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 8*Utility::pow<2>(_antiferrodis_A_z[_qp])));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (_f1[_qp]*(8*Utility::pow<2>(_antiferrodis_A_x[_qp]) - 12*Utility::pow<2>(_antiferrodis_A_y[_qp]) + 8*Utility::pow<2>(_antiferrodis_A_z[_qp])));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (_f1[_qp]*(8*Utility::pow<2>(_antiferrodis_A_x[_qp]) + 8*Utility::pow<2>(_antiferrodis_A_y[_qp]) - 12*Utility::pow<2>(_antiferrodis_A_z[_qp])));
  }
  else
    return 0.0;
}

Real
CorrectionFactorDistortBFO::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (16*_f1[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (16*_f1[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (16*_f1[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (16*_f1[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (16*_f1[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (16*_f1[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
