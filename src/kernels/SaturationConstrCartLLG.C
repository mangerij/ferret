/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the ter_Ms of the GNU General Public License _As published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET ple_Ase contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "SaturationConstrCartLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", SaturationConstrCartLLG);

template<>
InputParameters validParams<SaturationConstrCartLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual (penalty) contribution for saturation energy. If |m| = M_s then this energy should be zero.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("As", "As");
  params.addRequiredParam<Real>("Ms", "Ms");
  return params;
}

SaturationConstrCartLLG::SaturationConstrCartLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _alpha(getParam<Real>("alpha")),
  _g0(getParam<Real>("g0")),
  _As(getParam<Real>("As")),
  _Ms(getParam<Real>("Ms"))
{
}

Real
SaturationConstrCartLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (2*_alpha*_As*_g0*_mag_x[_qp]*(-1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))*(-_Ms + std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])));
  }
  else if (_component == 1)
  {
    return (2*_alpha*_As*_g0*_mag_y[_qp]*(-1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))*(-_Ms + std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])));
  }
  else if (_component == 2)
  {
    return (2*_alpha*_As*_g0*_mag_z[_qp]*(-1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))*(-_Ms + std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])));
  }
  else
    return 0.0;
}

Real
SaturationConstrCartLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (2*_alpha*_As*_g0*(Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*(-1 + 3*Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) - 
       _Ms*(2*Utility::pow<4>(_mag_x[_qp]) + 3*Utility::pow<2>(_mag_x[_qp])*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) + (-1 + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
  }
  else if (_component == 1)
  {
    return (2*_alpha*_As*_g0*(Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*(-1 + Utility::pow<2>(_mag_x[_qp]) + 3*Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) - 
       _Ms*(Utility::pow<4>(_mag_x[_qp]) + 2*Utility::pow<4>(_mag_y[_qp]) + (-1 + 3*Utility::pow<2>(_mag_y[_qp]))*Utility::pow<2>(_mag_z[_qp]) + Utility::pow<4>(_mag_z[_qp]) + Utility::pow<2>(_mag_x[_qp])*(-1 + 3*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
  }
  else if (_component == 2)
  {
    return (2*_alpha*_As*_g0*(Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*(-1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + 3*Utility::pow<2>(_mag_z[_qp])) - 
       _Ms*(Utility::pow<4>(_mag_x[_qp]) + Utility::pow<4>(_mag_y[_qp]) + 2*Utility::pow<4>(_mag_z[_qp]) + Utility::pow<2>(_mag_y[_qp])*(-1 + 3*Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(_mag_x[_qp])*(-1 + 2*Utility::pow<2>(_mag_y[_qp]) + 3*Utility::pow<2>(_mag_z[_qp]))))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
  }
  else
    return 0.0;
}

Real
SaturationConstrCartLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (2*_alpha*_As*_g0*_mag_x[_qp]*_mag_y[_qp]*(2*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) - _Ms*(1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
    }
    else if (jvar == _mag_z_var)
    {
      return (2*_alpha*_As*_g0*_mag_x[_qp]*_mag_z[_qp]*(2*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) - _Ms*(1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _mag_x_var)
    {
      return (2*_alpha*_As*_g0*_mag_x[_qp]*_mag_y[_qp]*(2*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) - _Ms*(1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
    }
    else if (jvar == _mag_z_var)
    {
      return (2*_alpha*_As*_g0*_mag_y[_qp]*_mag_z[_qp]*(2*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) - _Ms*(1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _mag_x_var)
    {
      return (2*_alpha*_As*_g0*_mag_x[_qp]*_mag_z[_qp]*(2*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) - _Ms*(1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
    }
    else if (jvar == _mag_y_var)
    {
      return (2*_alpha*_As*_g0*_mag_y[_qp]*_mag_z[_qp]*(2*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) - _Ms*(1 + Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha))*Utility::pow<3>(std::sqrt(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
