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

#ifndef TOTALENERGYALL_H
#define TOTALENERGYALL_H

//TODO: include the base header
#include "GeneralPostprocessor.h"

//Forward Declarations
class TotalEnergyAll;

template<>
InputParameters validParams<TotalEnergyAll>();

//TODO: change the base class!
class TotalEnergyAll : public GeneralPostprocessor
{
public:
  TotalEnergyAll(const InputParameters & parameters);
  virtual ~TotalEnergyAll();
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
protected:
  const PostprocessorValue & _Fbulk, & _Fwall, & _Felec, & _Fcoupled, & _Felastic, & _Fextelec;
};

#endif
