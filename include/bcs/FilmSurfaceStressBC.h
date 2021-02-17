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

#ifndef FILMSURFACESTRESSBC_H
#define FILMSURFACESTRESSBC_H

#include "IntegratedBC.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class FilmSurfaceStressBC;

template<>
InputParameters validParams<FilmSurfaceStressBC>();

class FilmSurfaceStressBC : public IntegratedBC
{
public:
  FilmSurfaceStressBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  std::string _base_name;
  const int _component;

  const MaterialProperty<RankTwoTensor> & _stress;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _C11;
  const MaterialProperty<Real> & _C12;
  const MaterialProperty<Real> & _C44;
  const MaterialProperty<Real> & _Q11;
  const MaterialProperty<Real> & _Q12;
  const MaterialProperty<Real> & _Q44;
};

#endif
