/**
 * @file   PiezoelectricApprox.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate piezoelectric constant 
 * d_{33} = P_z / \sigma_{zz}
 * for more information, see IEEE Trans. Ultra. Ferro. Freq. Contr. 51, 3, (2004)
 *
 */


#include "PiezoelectricApprox.h"

template<>
InputParameters validParams<PiezoelectricApprox>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

PiezoelectricApprox::PiezoelectricApprox(const InputParameters & parameters) :
  AuxKernel(parameters ),
   _stress( getMaterialProperty<RankTwoTensor>("stress")),
   _polar_z(coupledValue("polar_z"))
{
}

Real
PiezoelectricApprox::computeValue()
{
    return _polar_z[_qp] / _stress[_qp](2,2);
    
}


