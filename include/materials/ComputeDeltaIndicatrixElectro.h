/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#ifndef COMPUTEDELTAINDICATRIXELECTRO_H
#define COMPUTEDELTAINDICATRIXELECTRO_H

#include "RankThreeTensor.h"
#include "ComputeElectroopticTensorBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputeDeltaIndicatrixElectro defines an impermeability tensor material object with a given base name.
 */
class ComputeDeltaIndicatrixElectro : public ComputeElectroopticTensorBase
{
public:
  ComputeDeltaIndicatrixElectro(const InputParameters & parameters);

protected:
  virtual void ComputeDeltaIndicatrixElectro();
  const MaterialProperty<Real> & _permittivity;
  const MaterialProperty<Real> & _len_scale;
  const MaterialProperty<RankThreeTensor> & _electrooptic_tensor;
};

#endif //COMPUTEDELTAINDICATRIXELECTRO_H
