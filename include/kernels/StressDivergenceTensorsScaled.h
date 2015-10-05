/**
 * @file   StressDivergenceTensorsScaled.h
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15 2015
 * @brief  StressDivergenceTensorsScaled allows manually scaling of the StressDivergerceTensors. A propert scaling is crucial for numerical robustness.
 */
#ifndef STRESSDIVERGENCETENSORSSCALED_H
#define STRESSDIVERGENCETENSORSSCALED_H

#include "StressDivergenceTensors.h"
class StressDivergenceTensorsScaled;
template<>
InputParameters validParams<StressDivergenceTensorsScaled>();

class StressDivergenceTensorsScaled : public StressDivergenceTensors
{
public:
  StressDivergenceTensorsScaled(const InputParameters & parameters);
  ~StressDivergenceTensorsScaled(){}

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const Real _strain_scale;  //scaling the strains as in Ng et al
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm

};
#endif //STRESSDIVERGENCETENSORSSCALED_H
