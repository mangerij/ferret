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

#ifndef ENERGYRATEPOSTPROCESSOR_H
#define ENERGYRATEPOSTPROCESSOR_H

#include "GeneralPostprocessor.h"

class EnergyRatePostprocessor;

template<>
InputParameters validParams<EnergyRatePostprocessor>();

/**
 * This postprocessor displays the change in the postprocessor between
 * adjacent timesteps divided by the time step
 */

 class EnergyRatePostprocessor : public GeneralPostprocessor
 {
 public:
   EnergyRatePostprocessor(const InputParameters & parameters);
   virtual void initialize();
   virtual void execute();
   virtual Real getValue();
 protected:
   const PostprocessorValue & _postprocessor, & _postprocessor_old, & _dt, & _dt_old;
 };

#endif /* ENERGYRATEPOSTPROCESSOR_H */
