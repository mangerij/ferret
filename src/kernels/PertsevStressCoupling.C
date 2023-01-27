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

#include "PertsevStressCoupling.h"

class PertsevStressCoupling;

registerMooseObject("FerretApp", PertsevStressCoupling);

InputParameters PertsevStressCoupling::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to a stress-formulation of the electrostrictive coupling.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("Q11", "the 11 (Voight) component of electrostrictive coupling coefficient");
  params.addRequiredParam<Real>("Q12", "the 12 (Voight) component of electrostrictive coupling coefficient");
  params.addRequiredParam<Real>("Q44", "the 44 (Voight) component of electrostrictive coupling coefficient");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

PertsevStressCoupling::PertsevStressCoupling(const InputParameters & parameters)
  :Kernel(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
   _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _Q11(getParam<Real>("Q11")),
   _Q12(getParam<Real>("Q12")),
   _Q44(getParam<Real>("Q44")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
PertsevStressCoupling::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (-2*_polar_x[_qp]*_Q11*_stress[_qp](0,0) - _Q44*(_polar_y[_qp]*_stress[_qp](0,1) + _polar_z[_qp]*_stress[_qp](0,2)) - _Q12*(2*_polar_x[_qp]*_stress[_qp](1,1) + 2*_polar_x[_qp]*_stress[_qp](2,2)));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (-2*_polar_y[_qp]*_Q11*_stress[_qp](1,1) - _Q44*(_polar_x[_qp]*_stress[_qp](0,1) + _polar_z[_qp]*_stress[_qp](1,2)) - _Q12*(2*_polar_y[_qp]*_stress[_qp](0,0) + 2*_polar_y[_qp]*_stress[_qp](2,2)));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (-(_Q12*(2*_polar_z[_qp]*_stress[_qp](0,0) + 2*_polar_z[_qp]*_stress[_qp](1,1))) - _Q44*(_polar_x[_qp]*_stress[_qp](0,2) + _polar_y[_qp]*_stress[_qp](1,2)) - 2*_polar_z[_qp]*_Q11*_stress[_qp](2,2));
  }
  else
    return 0.0;
}

Real
PertsevStressCoupling::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*_Q11*_stress[_qp](0,0) - _Q12*(2*_stress[_qp](1,1) + 2*_stress[_qp](2,2)));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*_Q11*_stress[_qp](1,1) - _Q12*(2*_stress[_qp](0,0) + 2*_stress[_qp](2,2)));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-(_Q12*(2*_stress[_qp](0,0) + 2*_stress[_qp](1,1))) - 2*_Q11*_stress[_qp](2,2));
  }
  else
    return 0.0;
}

Real
PertsevStressCoupling::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-(_Q44*_stress[_qp](0,1)));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-(_Q44*_stress[_qp](0,2)));
    }
    //else if (jvar == _disp_x_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    //else if (jvar == _disp_y_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    //else if (jvar == _disp_z_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-(_Q44*_stress[_qp](0,1)));
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-(_Q44*_stress[_qp](1,2)));
    }
    //else if (jvar == _disp_x_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    //else if (jvar == _disp_y_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    //else if (jvar == _disp_z_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-(_Q44*_stress[_qp](0,2)));
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-(_Q44*_stress[_qp](1,2)));
    }
    //else if (jvar == _disp_x_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    //else if (jvar == _disp_y_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    //else if (jvar == _disp_z_var)
    //{
    //  return _test[_i][_qp] * ();
    //}
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
