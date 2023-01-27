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

#ifndef COMPUTEPIEZOTENSOR_H
#define COMPUTEPIEZOTENSOR_H

#include "RankThreeTensor.h"
#include "RankFourTensor.h"
#include "PiezostrictiveTensorTools.h"
#include "ComputeRotatedPiezoTensorBase.h"

/**
 * ComputeElectroopticTensor defines a linear electrooptic tensor material object with a given base name.
 */
class ComputePiezoTensor : public ComputeRotatedPiezoTensorBase
{
public:
  ComputePiezoTensor(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpPiezoTensor();

  /// Individual material information
  bool _compute_piezostrictive_coeff;
  RankThreeTensor _dkpq;
  RankThreeTensor _dpqkT;
  RankThreeTensor _Dkij;
  RankThreeTensor _DijkT;
  RankFourTensor _Cijkl;
};

#endif //COMPUTEPIEZOTENSOR_H
