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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef CONVERSEPIEZOELECTRICSTRAIN_H
#define CONVERSEPIEZOELECTRICSTRAIN_H

#include "Kernel.h"
#include "ComputePiezostrictiveTensor.h"
#include "Material.h"

//Forward Declarations
class ConversePiezoelectricStrain;

template<>
InputParameters validParams<ConversePiezoelectricStrain>();

class ConversePiezoelectricStrain: public Kernel
{
public:
  ConversePiezoelectricStrain(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const MaterialProperty<RankThreeTensor> & _piezostrictive_tensor;
  const unsigned int _component;
  const unsigned int _potential_E_int_var;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm

};
#endif //CONVERSEPIEZOELECTRICSTRAIN_H
