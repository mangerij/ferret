/**
 * @file   ComputeIndicatrixBase.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate the indicatrix
 *
 */

#include "ComputeIndicatrixBase.h"

template<>
InputParameters validParams<ComputeIndicatrixBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");

  return params;
}

ComputeIndicatrixBase::ComputeIndicatrixBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _indicatrix_name(_base_name + "indicatrix"),
   _indicatrix(declareProperty<RealVectorValue>(_indicatrix_name))
{
}

void
ComputeIndicatrixBase::computeQpProperties()
{
  computeQpIndicatrix();
}
