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

#ifndef LANDAUMATERIALBULKFOURTHDERIVATIVE_H
#define LANDAUMATERIALBULKFOURTHDERIVATIVE_H

#include "Kernel.h"
#include "ComputeRankTwoLandauTensor.h"
#include "ComputeRankFourLandauTensor.h"

class LandauMaterialBulkFourthDerivative;

template<>
InputParameters validParams<LandauMaterialBulkFourthDerivative>();

class LandauMaterialBulkFourthDerivative: public Kernel
{
public:

  LandauMaterialBulkFourthDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<RankTwoTensor> & _rank_two_landau_tensor;
  const MaterialProperty<RankFourTensor> & _rank_four_landau_tensor;
  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _len_scale;
};
#endif //LANDAUMATERIALBULKFOURTHDERIVATIVE_H