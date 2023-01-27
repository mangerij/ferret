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

#ifndef SURFACEMECHANICSBC_H
#define SURFACEMECHANICSBC_H

#include "IntegratedBC.h"

template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
template <typename>
class RankFourTensorTempl;
typedef RankFourTensorTempl<Real> RankFourTensor;

class SurfaceMechanicsBC : public IntegratedBC
{
public:

  SurfaceMechanicsBC(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

  virtual RankTwoTensor computeQpProjection();

  virtual RankTwoTensor computeQpRankTwoRotation(unsigned int tensor_component);

  virtual RankFourTensor computeQpRankFourRotation(unsigned int tensor_component);

  const unsigned int _component;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;

  std::vector<Real> _Csijkl_vector; //surfaceFillFromInputVector method only takes std::vector<Real>. Might need to change..

  RealVectorValue _S_k_vector;

  Real _taus;

  RealVectorValue _surface_euler_angles;

};
#endif // SURFACEMECHANICSBC_H
