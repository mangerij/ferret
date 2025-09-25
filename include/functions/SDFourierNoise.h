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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

// adapted from FourierNoise.C [credit to D. Schwen (INL)]

#pragma once

#include "Function.h"

class FEProblemBase;

/**
 * Generate noise using random fourier series coefficients
 */
class SDFourierNoise : public Function
{
public:
  SDFourierNoise(const InputParameters & parameters);

  static InputParameters validParams();

  virtual Real value(Real, const Point & p) const override;

protected:
  struct SeriesItem
  {
    /// k-vector
    RealVectorValue k;
    /// sin coefficient
    Real s;
    /// cos coefficient
    Real c;
  };

  /// selected lower lengthscale for the noise cut-off
  const Real _lambda;

  /// Fourier series terms
  std::vector<SeriesItem> _series;

  /// amplitude factor
  Real _scale;

  /// actual range
  Real _range;

  /// actual mid
  Real _mid;

  /// FEProblem pointer for obtaining the current mesh
  FEProblemBase & _fe_problem;
};
