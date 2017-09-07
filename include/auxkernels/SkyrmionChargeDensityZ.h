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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

*/

#ifndef SKYRMIONCHARGEDENSITYZ_H
#define SKYRMIONCHARGEDENSITYZ_H

#include "AuxKernel.h"

//Forward declarations
class SkyrmionChargeDensityZ;

template<>
InputParameters validParams<SkyrmionChargeDensityZ>();

class SkyrmionChargeDensityZ : public AuxKernel
{
public:
  SkyrmionChargeDensityZ(const InputParameters & parameters);

  virtual ~SkyrmionChargeDensityZ() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue & _polar_x_norm;
  const VariableValue & _polar_y_norm;
  const VariableValue & _polar_z_norm;
  const VariableGradient & _polar_x_norm_grad;
  const VariableGradient & _polar_y_norm_grad;
  const VariableGradient & _polar_z_norm_grad;
};

#endif
