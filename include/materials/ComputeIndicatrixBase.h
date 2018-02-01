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

#ifndef COMPUTEINDICATRIXBASE_H
#define COMPUTEINDICATRIXBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

class ComputeIndicatrixBase;

template<>
InputParameters validParams<ComputeIndicatrixBase>();

/**
 * ComputeTensorBase the base class for computing photostrictive tensors
 */
class ComputeIndicatrixBase : public Material
{
public:
  ComputeIndicatrixBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpIndicatrix() = 0;

  std::string _base_name;
  std::string _indicatrix_name;

  MaterialProperty<RankTwoTensor> & _indicatrix;

};

#endif //COMPUTEINDICATRIXBASE_H
