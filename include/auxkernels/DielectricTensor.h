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

#ifndef DIELECTRICTENSOR_H
#define DIELECTRICTENSOR_H

#include "AuxKernel.h"
#include "ComputeElectrostrictiveTensor.h"

template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;

class DielectricTensor : public AuxKernel
{
public:
  DielectricTensor(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeValue();
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123;
  const unsigned int _first_deriv, _second_deriv;
  const Real _len_scale;
};

#endif // DIELECTRICTENSOR_H
