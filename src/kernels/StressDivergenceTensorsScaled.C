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
  params.addParam<Real>("energy_scale",1.0,"energy scale");
  params.addParam<Real>("disp_scale", 1.0, "scale for the displacement variables");
  params.addParam<Real>("elastic_equation_scale",1.0,"other scaling in the residual equation for elasticity variables");
  return params;
}

StressDivergenceTensorsScaled::StressDivergenceTensorsScaled(const std::string& name, InputParameters parameters)
  :StressDivergenceTensors(name, parameters),
   _len_scale(getParam<Real>("len_scale")),
   _energy_scale(getParam<Real>("energy_scale")),
   _disp_scale(getParam<Real>("disp_scale")),
   _elastic_equation_scale(getParam<Real>("elastic_equation_scale"))
{

}
Real
StressDivergenceTensorsScaled::computeQpResidual()
{
  return _elastic_equation_scale*pow(_disp_scale,2)*_energy_scale*_len_scale*StressDivergenceTensors::computeQpResidual();
}

Real
StressDivergenceTensorsScaled::computeQpJacobian()
{
  return  _elastic_equation_scale*pow(_disp_scale,2)*_energy_scale*_len_scale*StressDivergenceTensors::computeQpJacobian();
}

Real
StressDivergenceTensorsScaled::computeQpOffDiagJacobian(unsigned int jvar)
{
  return  _elastic_equation_scale*pow(_disp_scale,2)*_energy_scale*_len_scale*StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}
