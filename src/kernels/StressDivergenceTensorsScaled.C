/**
 * @file   StressDivergenceTensorsScaled.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Nov 25 13:10:24 2013
 * @brief  StressDivergenceTensorsScaled allows manually scaling of the StressDivergerceTensors. A propert scaling is crucial for numerical robustness.
 */

#include "StressDivergenceTensorsScaled.h"

template<>
InputParameters validParams<StressDivergenceTensorsScaled>()
{
  InputParameters params = validParams<StressDivergenceTensors>();
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

StressDivergenceTensorsScaled::StressDivergenceTensorsScaled(const std::string& name, InputParameters parameters)
  :StressDivergenceTensors(name, parameters),
   _len_scale(getParam<Real>("len_scale"))
{

}
Real
StressDivergenceTensorsScaled::computeQpResidual()
{
  return _len_scale*StressDivergenceTensors::computeQpResidual();
}

Real
StressDivergenceTensorsScaled::computeQpJacobian()
{
  return  _len_scale*StressDivergenceTensors::computeQpJacobian();
}

Real
StressDivergenceTensorsScaled::computeQpOffDiagJacobian(unsigned int jvar)
{
  return  _len_scale*StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}
