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

#include "DomainVariantPopulation.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", DomainVariantPopulation);

InputParameters DomainVariantPopulation::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates the fraction of volume of a given domain population (only works in tetragonal phase at the moment)");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<PostprocessorName>("vol", "The name of the volume used for volume normalization");
  return params;
}

DomainVariantPopulation::DomainVariantPopulation(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _component(getParam<unsigned int>("component")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _vol(getPostprocessorValue("vol"))
{
}

Real
DomainVariantPopulation::computeQpIntegral()
{
if (_component == 0)
  {
    return (1.0/_vol)*std::abs(_polar_x[_qp])/std::sqrt(_polar_x[_qp]*_polar_x[_qp]+_polar_y[_qp]*_polar_y[_qp]+_polar_z[_qp]*_polar_z[_qp]);
  }
  else if (_component == 1)
  {
    return (1.0/_vol)*std::abs(_polar_y[_qp])/std::sqrt(_polar_x[_qp]*_polar_x[_qp]+_polar_y[_qp]*_polar_y[_qp]+_polar_z[_qp]*_polar_z[_qp]);
  }
  else if (_component == 2)
  {
    return (1.0/_vol)*std::abs(_polar_z[_qp])/std::sqrt(_polar_x[_qp]*_polar_x[_qp]+_polar_y[_qp]*_polar_y[_qp]+_polar_z[_qp]*_polar_z[_qp]);
  }
  else
    return 0.0;
}
