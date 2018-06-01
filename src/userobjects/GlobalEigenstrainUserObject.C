/**
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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "GlobalEigenstrainUserObject.h"
#include "GlobalStrainUserObject.h"

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", GlobalEigenstrainUserObject);

template <>
InputParameters
validParams<GlobalEigenstrainUserObject>()
{
  InputParameters params = validParams<GlobalStrainUserObject>();
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the afd vector field");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("C11", "the 11 component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C12", "the 12 component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C44", "the 44 component of elastic stiffness tensor");
  params.addRequiredParam<Real>("Q11", "the 11 component of electrostrictive coupling tensor");
  params.addRequiredParam<Real>("Q12", "the 12 component of electrostrictive coupling tensor");
  params.addRequiredParam<Real>("Q44", "the 44 component of electrostrictive coupling tensor");
  params.addRequiredParam<Real>("R11", "the 11 component of rotostrictive coupling tensor");
  params.addRequiredParam<Real>("R12", "the 12 component of rotostrictive coupling tensor");
  params.addRequiredParam<Real>("R44", "the 44 component of rotostrictive coupling tensor");
  return params;
}

GlobalEigenstrainUserObject::GlobalEigenstrainUserObject(const InputParameters & parameters)
  : GlobalStrainUserObject(parameters),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _C11(getParam<Real>("C11")),
   _C12(getParam<Real>("C12")),
   _C44(getParam<Real>("C44")),
   _Q11(getParam<Real>("Q11")),
   _Q12(getParam<Real>("Q12")),
   _Q44(getParam<Real>("Q44")),
   _R11(getParam<Real>("R11")),
   _R12(getParam<Real>("R12")),
   _R44(getParam<Real>("R44"))
{
}

void
GlobalStrainUserObject::initialize()
{
}

void
GlobalEigenstrainUserObject::computeAdditionalStress()
{ 
  for (unsigned int _qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    _applied_stress_tensor(0,0) += _C11*Utility::pow<2>(_polar_x[_qp])*_Q11 + _C12*Utility::pow<2>(_polar_y[_qp])*_Q11 + _C12*Utility::pow<2>(_polar_z[_qp])*_Q11 
			        + 2.0*_C12*Utility::pow<2>(_polar_x[_qp]) * _Q12 + _C11*Utility::pow<2>(_polar_y[_qp])*_Q12 + _C12*Utility::pow<2>(_polar_y[_qp])*_Q12 
                                + _C11*Utility::pow<2>(_polar_z[_qp])*_Q12 + _C12*Utility::pow<2>(_polar_z[_qp])*_Q12 + _C11*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R11 
                                + _C12*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R11 + _C12*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R11 + 2.0*_C12*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12 
                                + _C11*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12 + _C12*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12 + _C11*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12 
                                + _C12*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12;

    _applied_stress_tensor(1,1) += _C12*Utility::pow<2>(_polar_x[_qp])*_Q11 + _C11*Utility::pow<2>(_polar_y[_qp])*_Q11 + _C12*Utility::pow<2>(_polar_z[_qp])*_Q11 
                                + _C11*Utility::pow<2>(_polar_x[_qp])*_Q12 + _C12*Utility::pow<2>(_polar_x[_qp])*_Q12 + 2.0*_C12*Utility::pow<2>(_polar_y[_qp])*_Q12 
                                + _C11*Utility::pow<2>(_polar_z[_qp])*_Q12 + _C12*Utility::pow<2>(_polar_z[_qp])*_Q12 + _C12*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R11 
                                + _C11*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R11 + _C12*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R11 + _C11*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12 
                                + _C12*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12 + 2.0*_C12*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12 + _C11*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12 
                                + _C12*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12;

    _applied_stress_tensor(2,2) += _C12*Utility::pow<2>(_polar_x[_qp])*_Q11 + _C12*Utility::pow<2>(_polar_y[_qp])*_Q11 + _C11*Utility::pow<2>(_polar_z[_qp])*_Q11 
                                + _C11*Utility::pow<2>(_polar_x[_qp])*_Q12 + _C12*Utility::pow<2>(_polar_x[_qp])*_Q12 + _C11*Utility::pow<2>(_polar_y[_qp])*_Q12 
                                + _C12*Utility::pow<2>(_polar_y[_qp])*_Q12 + 2.0*_C12*Utility::pow<2>(_polar_z[_qp])*_Q12 + _C12*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R11 
                                + _C12*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R11 + _C11*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R11 + _C11*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12 
                                + _C12*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12 + _C11*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12 + _C12*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12 
                                + 2.0*_C12*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12;

    _applied_stress_tensor(0,1) += 4.0*_C44*_polar_x[_qp]*_polar_y[_qp]*_Q44 + 4.0*_C44*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_R44;

    _applied_stress_tensor(1,2) += 4.0*_C44*_polar_y[_qp]*_polar_z[_qp]*_Q44 + 4.0*_C44*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_R44;

    _applied_stress_tensor(0,2) += 4.0*_C44*_polar_x[_qp]*_polar_z[_qp]*_Q44 + 4.0*_C44*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_R44;
  }
  _applied_stress_tensor(1,0) = _applied_stress_tensor(0,1);
  _applied_stress_tensor(2,1) = _applied_stress_tensor(1,2);
  _applied_stress_tensor(2,0) = _applied_stress_tensor(0,2);
}


void
GlobalStrainUserObject::execute()
{
}

void
GlobalStrainUserObject::threadJoin(const UserObject & uo)
{
}

void
GlobalStrainUserObject::finalize()
{
}


