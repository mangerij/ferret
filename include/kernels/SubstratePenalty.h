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

#ifndef SUBSTRATEPENALTY_H
#define SUBSTRATEPENALTY_H

#include "Kernel.h"

//Forward Declarations
class SubstratePenalty;

template<>
InputParameters validParams<SubstratePenalty>();

class SubstratePenalty: public Kernel
{
public:

  SubstratePenalty(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const unsigned int _component;
  const PostprocessorValue & _aveexx;
  const PostprocessorValue & _aveeyy;
  const VariableGradient & _u_x_grad;
  const VariableGradient & _u_y_grad;
  const VariableGradient & _u_z_grad;
  const MaterialProperty<Real> & _K;
  const MaterialProperty<Real> & _e0xx;
  const MaterialProperty<Real> & _e0yy;

};
#endif //SUBSTRATEPENALTY_H
