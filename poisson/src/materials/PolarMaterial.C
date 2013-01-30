// Original class author: S. Gu

#include "PolarMaterial.h"

/**
 * This class defines the polarization vector field for the Poisson problem
 **/

template<>
InputParameters validParams<PolarMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<RealVectorValue>("P", "Polarization");  
  return params;
}

PolarMaterial::PolarMaterial(const std::string & name, InputParameters parameters)
    : Material(name, parameters),
      _p(declareProperty<RealVectorValue>("Polarization")),
      _p_in(getParam<RealVectorValue>("P"))
      
{
  //std::cout<<"("<<_p_in(0)<<","<<_p_in(1)<<","<<_p_in(2)<<")";
  //std::cout<<"\n"; 
}

void
PolarMaterial::computeQpProperties()
{
  _p[_qp]=_p_in;
}
