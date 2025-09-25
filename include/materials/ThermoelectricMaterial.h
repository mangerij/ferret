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

#pragma once

#include "Material.h"
#include "MooseTypes.h"

// Forward Declarations
class Function;

/**
 * Simple material with constant properties.
 */
template <bool is_ad>
class ThermoelectricMaterialTempl : public Material
{
public:
  static InputParameters validParams();

  ThermoelectricMaterialTempl(const InputParameters & parameters);

protected:
  virtual void computeProperties();

  const bool _has_temp;
  const VariableValue & _T;
  const ADVariableValue & _ad_T;

  const Real _my_thC;
  const Real _my_ecC;
  const Real _my_sbC;

  GenericMaterialProperty<Real, is_ad> & _thC;
  MaterialProperty<Real> & _thC_dT;
  const Function * const _thC_temperature_function;

  GenericMaterialProperty<Real, is_ad> & _ecC;
  MaterialProperty<Real> & _ecC_dT;
  const Function * const _ecC_temperature_function;

  GenericMaterialProperty<Real, is_ad> & _sbC;
  MaterialProperty<Real> & _sbC_dT;
  const Function * const _sbC_temperature_function;

private:
  void setDerivatives(GenericReal<is_ad> & prop, Real dprop_dT, const ADReal & ad_T);
};

typedef ThermoelectricMaterialTempl<false> ThermoelectricMaterial;
typedef ThermoelectricMaterialTempl<true> ADThermoelectricMaterial;
