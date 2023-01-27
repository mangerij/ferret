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

#include "ConvertField.h"

registerMooseObject("FerretApp", ConvertField);

InputParameters ConvertField::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("This file is useful for quickly converting from reduced to real units in BFO based potentials.");
  params.addRequiredParam<unsigned int>("conv_type","type of conversion");
  params.addRequiredParam<unsigned int>("component","component of field to convert");
  params.addRequiredCoupledVar("polar_x_red", "The reduced field");
  params.addRequiredCoupledVar("polar_y_red", "The reduced field");
  params.addCoupledVar("polar_z_red", 0.0, "The reduced field");
  params.addRequiredCoupledVar("antiferrodis_A_x_red", "The reduced field");
  params.addRequiredCoupledVar("antiferrodis_A_y_red", "The reduced field");
  params.addCoupledVar("antiferrodis_A_z_red", 0.0, "The reduced field");
  params.addRequiredCoupledVar("disp_x_red", "The reduced field");
  params.addRequiredCoupledVar("disp_y_red", "The reduced field");
  params.addCoupledVar("disp_z_red", 0.0, "The reduced field");
  params.addRequiredCoupledVar("strain_xx_red", "The reduced field");
  params.addRequiredCoupledVar("strain_yy_red", "The reduced field");
  params.addCoupledVar("strain_zz_red", 0.0, "The reduced field");
  params.addRequiredCoupledVar("strain_xy_red", "The reduced field");
  params.addCoupledVar("strain_xz_red", 0.0, "The reduced field");
  params.addCoupledVar("strain_yz_red", 0.0, "The reduced field");
  params.addRequiredParam<Real>("Ps","spontaneous polarization");
  params.addRequiredParam<Real>("As","spontaneous tilt");
  params.addRequiredParam<Real>("uscale","spontaneous displacements scaling");
  params.addRequiredParam<Real>("es_norm","spontaneous normal strain");
  params.addRequiredParam<Real>("es_shear","spontaneous shear strain");
  return params;
}

ConvertField::ConvertField(const InputParameters & parameters) :
  AuxKernel(parameters),
   _conv_type(getParam<unsigned int>("conv_type")),
   _component(getParam<unsigned int>("component")),
   _polar_x_red(coupledValue("polar_x_red")),
   _polar_y_red(coupledValue("polar_y_red")),
   _polar_z_red(coupledValue("polar_z_red")),
   _antiferrodis_A_x_red(coupledValue("antiferrodis_A_x_red")),
   _antiferrodis_A_y_red(coupledValue("antiferrodis_A_y_red")),
   _antiferrodis_A_z_red(coupledValue("antiferrodis_A_z_red")),
   _disp_x_red(coupledValue("disp_x_red")),
   _disp_y_red(coupledValue("disp_y_red")),
   _disp_z_red(coupledValue("disp_z_red")),
   _strain_xx_red(coupledValue("strain_xx_red")),
   _strain_yy_red(coupledValue("strain_yy_red")),
   _strain_zz_red(coupledValue("strain_zz_red")),
   _strain_xy_red(coupledValue("strain_xz_red")),
   _strain_xz_red(coupledValue("strain_xz_red")),
   _strain_yz_red(coupledValue("strain_yz_red")),
   _Ps(getParam<Real>("Ps")),
   _As(getParam<Real>("As")),
   _uscale(getParam<Real>("uscale")),
   _es_norm(getParam<Real>("es_norm")),
   _es_shear(getParam<Real>("es_shear"))
{}

Real
ConvertField::computeValue()
{
  if (_conv_type == 0)
  {
    if (_component == 0)
    {
      return _uscale*_disp_x_red[_qp];
    }
    else if (_component == 1)
    {
      return _uscale*_disp_y_red[_qp];
    }
    else if (_component == 2)
    {
      return _uscale*_disp_z_red[_qp];
    }
    else
      return 0.0;
  }
  else if (_conv_type == 1)
  {
    if (_component == 0)
    {
      return _Ps*_polar_x_red[_qp];
    }
    else if (_component == 1)
    {
      return _Ps*_polar_y_red[_qp];
    }
    else if (_component == 2)
    {
      return _Ps*_polar_z_red[_qp];
    }
    else
      return 0.0;
  }
  else if (_conv_type == 2)
  {
    if (_component == 0)
    {
      return _As*_antiferrodis_A_x_red[_qp];
    }
    else if (_component == 1)
    {
      return _As*_antiferrodis_A_y_red[_qp];
    }
    else if (_component == 2)
    {
      return _As*_antiferrodis_A_z_red[_qp];
    }
    else
      return 0.0;
  }
  else if (_conv_type == 3)
  {
    if (_component == 0)
    {
      return _es_norm*_strain_xx_red[_qp];
    }
    else if (_component == 1)
    {
      return _es_norm*_strain_yy_red[_qp];
    }
    else if (_component == 2)
    {
      return _es_norm*_strain_zz_red[_qp];
    }
    else
      return 0.0;
  }
  else if (_conv_type == 4)
  {
    if (_component == 0)
    {
      return _es_shear*_strain_xy_red[_qp];
    }
    else if (_component == 1)
    {
      return _es_shear*_strain_xz_red[_qp];
    }
    else if (_component == 2)
    {
      return _es_shear*_strain_yz_red[_qp];
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
