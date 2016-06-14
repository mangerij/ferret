/**
 * @file   StressDivergenceTensorsScaled.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Mon Nov 25 13:10:24 2013
 * @brief  StressDivergenceTensorsScaled allows manually scaling of the StressDivergerceTensors.
 */

#include "StressDivergenceTensorsScaled.h"

template<>
InputParameters validParams<StressDivergenceTensorsScaled>()
{
  InputParameters params = validParams<StressDivergenceTensors>();
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

StressDivergenceTensorsScaled::StressDivergenceTensorsScaled(const InputParameters & parameters)
  :StressDivergenceTensors(parameters),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
StressDivergenceTensorsScaled::computeQpResidual()
{
  return std::pow(_len_scale, 2.0) * StressDivergenceTensors::computeQpResidual();
}

Real
StressDivergenceTensorsScaled::computeQpJacobian()
{
  return  std::pow(_len_scale, 2.0) * StressDivergenceTensors::computeQpJacobian();
}

Real
StressDivergenceTensorsScaled::computeQpOffDiagJacobian(unsigned int jvar)
{
  return  std::pow(_len_scale, 2.0) * StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}
