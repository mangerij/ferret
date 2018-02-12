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

#ifndef COMPUTERANKFOURLANDAUTENSORBASE_H
#define COMPUTERANKFOURLANDAUTENSORBASE_H

#include "Material.h"
#include "RankFourTensor.h"

class ComputeRankFourLandauTensorBase;

template<>
InputParameters validParams<ComputeRankFourLandauTensorBase>();

/**
 * ComputePiezoTensorBase the base class for computing photostrictive tensors
 */
class ComputeRankFourLandauTensorBase : public Material
{
public:
  ComputeRankFourLandauTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpRankFourLandauTensor() = 0;

  std::string _base_name;
  std::string _rank_four_landau_tensor_name;
  MaterialProperty<RankFourTensor> & _rank_four_landau_tensor;
};

#endif //COMPUTERANKFOURLANDAUTENSORBASE_H
