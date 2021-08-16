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

#include "CorrectionFactorCoupledDistortBFO.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CorrectionFactorCoupledDistortBFO);

template<>
InputParameters validParams<CorrectionFactorCoupledDistortBFO>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a correction factor for the local free energy");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt");
  return params;
}

CorrectionFactorCoupledDistortBFO::CorrectionFactorCoupledDistortBFO(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _c0(getMaterialProperty<Real>("c0"))
{
}

Real
CorrectionFactorCoupledDistortBFO::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (0.0);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (0.0);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (0.0);
  }
  else
    return 0.0;
}

Real
CorrectionFactorCoupledDistortBFO::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
  }
  else
    return 0.0;
}

Real
CorrectionFactorCoupledDistortBFO::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
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
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.0);
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
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.0);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.0);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
