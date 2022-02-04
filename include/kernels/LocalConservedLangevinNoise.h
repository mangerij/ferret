//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LocalLangevinNoise.h"
#include "LocalConservedNoiseBase.h"

class LocalConservedLangevinNoise : public LocalLangevinNoise
{
public:
  static InputParameters validParams();

  LocalConservedLangevinNoise(const InputParameters & parameters);

protected:
  virtual void residualSetup(){};
  virtual Real computeQpResidual();

private:
  const LocalConservedNoiseInterface & _noise;

  const VariableValue & _T_bath;
};
