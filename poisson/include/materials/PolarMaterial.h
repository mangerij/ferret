// Original class author: S. Gu

#ifndef POLARMATERIAL_H
#define POLARMATERIAL_H

#include "Material.h"

/**
 * This class defines the polarization vector field for the Poisson problem
 **/

class PolarMaterial;

template<>
InputParameters validParams<PolarMaterial>();

class PolarMaterial : public Material
{
public:
  PolarMaterial(const std:: string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();
  MaterialProperty<RealVectorValue >&  _p; 
  //vector to get input values
  RealVectorValue _p_in;
};

#endif //POLARMATERIAL_H
