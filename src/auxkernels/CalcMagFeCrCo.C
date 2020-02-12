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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "CalcMagFeCrCo.h"

registerMooseObject("FerretApp", CalcMagFeCrCo);

template<>

InputParameters validParams<CalcMagFeCrCo>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates concentration dependent magnetization field");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addCoupledVar("c1", "The first concentration variable");
  params.addCoupledVar("c2", "The second concentration variable");
  params.addCoupledVar("c3", "The third concentration variable");
  params.addRequiredParam<Real>("bohrM", "bohrM");
  params.addRequiredParam<Real>("T", "T");
  return params;
}


CalcMagFeCrCo::CalcMagFeCrCo(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
  _c1(coupledValue("c1")),
  _c2(coupledValue("c2")),
  _c3(coupledValue("c3")),
  _bohrM(getParam<Real>("bohrM")),
  _T(getParam<Real>("T"))
{
}

Real
CalcMagFeCrCo::computeValue()
{
  if (_component == 0)
  {
   Real tau = _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + 1450*_c3[_qp] + (1650 + 550*(_c2[_qp]-_c1[_qp]))*_c1[_qp]*_c2[_qp] + 590*_c1[_qp]*_c3[_qp]);
   if (tau > 0.9)
   {
     return _bohrM*(2.22*_c1[_qp] - 0.01*_c2[_qp] + 1.35*_c3[_qp] - 0.85*_c1[_qp]*_c2[_qp] + (2.4127 + 0.2418*(_c3[_qp]-_c1[_qp]))*_c1[_qp]*_c3[_qp])*std::pow(2.0,-(2.0+10.0*(tau-1.0)));
   }
   else if (tau <= 0.9)
   {
     return _bohrM*(2.22*_c1[_qp] - 0.01*_c2[_qp] + 1.35*_c3[_qp] - 0.85*_c1[_qp]*_c2[_qp] + (2.4127 + 0.2418*(_c3[_qp]-_c1[_qp]))*_c1[_qp]*_c3[_qp])*(1.0 - (1.0/7.0)*(5.0*std::pow(tau,4.0)+2.0*std::pow(tau,20.0)));
   }
   else
     return 0.0;
  }
  else if (_component == 1)
  {
   Real tau = _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + 1450*_c3[_qp] + (1650 + 550*(_c2[_qp]-_c1[_qp]))*_c1[_qp]*_c2[_qp] + 590*_c1[_qp]*_c3[_qp]);
   if (tau > 0.9)
   {
     return _bohrM*(2.22*_c1[_qp] - 0.01*_c2[_qp] + 1.35*_c3[_qp] - 0.85*_c1[_qp]*_c2[_qp] + (2.4127 + 0.2418*(_c3[_qp]-_c1[_qp]))*_c1[_qp]*_c3[_qp])*std::pow(2.0,-(2.0+10.0*(tau-1.0)));
   }
   else if (tau <= 0.9)
   {
     return _bohrM*(2.22*_c1[_qp] - 0.01*_c2[_qp] + 1.35*_c3[_qp] - 0.85*_c1[_qp]*_c2[_qp] + (2.4127 + 0.2418*(_c3[_qp]-_c1[_qp]))*_c1[_qp]*_c3[_qp])*(1.0 - (1.0/7.0)*(5.0*std::pow(tau,4.0)+2.0*std::pow(tau,20.0)));
   }
   else
     return 0.0;
  }
  else if (_component == 2)
  {
   Real tau = _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + 1450*_c3[_qp] + (1650 + 550*(_c2[_qp]-_c1[_qp]))*_c1[_qp]*_c2[_qp] + 590*_c1[_qp]*_c3[_qp]);
   if (tau > 0.9)
   {
     return _bohrM*(2.22*_c1[_qp] - 0.01*_c2[_qp] + 1.35*_c3[_qp] - 0.85*_c1[_qp]*_c2[_qp] + (2.4127 + 0.2418*(_c3[_qp]-_c1[_qp]))*_c1[_qp]*_c3[_qp])*std::pow(2.0,-(2.0+10.0*(tau-1.0)));
   }
   else if (tau <= 0.9)
   {
     return _bohrM*(2.22*_c1[_qp] - 0.01*_c2[_qp] + 1.35*_c3[_qp] - 0.85*_c1[_qp]*_c2[_qp] + (2.4127 + 0.2418*(_c3[_qp]-_c1[_qp]))*_c1[_qp]*_c3[_qp])*(1.0 - (1.0/7.0)*(5.0*std::pow(tau,4.0)+2.0*std::pow(tau,20.0)));
   }
   else
     return 0.0;
  }
  else
    return 0.0;
}
