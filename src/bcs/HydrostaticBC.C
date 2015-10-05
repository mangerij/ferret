/****************************************************************/
/* Hydrostatic BC:                                              */
/*     This BC is intended only for testing purpose.            */
/****************************************************************/


#include "HydrostaticBC.h"

template<>
InputParameters validParams<HydrostaticBC>()
{
    InputParameters params = validParams<IntegratedBC>();
  // Here we are adding a parameter that will be extracted from the input file by the Parser
  params.addRequiredParam<Real>("pressure","Specify the hydrostatic pressure");
  params.addRequiredParam<int>("component","Component of displacement for BC");

  return params;
}

HydrostaticBC::HydrostaticBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _pressure(getParam<Real>("pressure")),
  _component(getParam<int>("component"))
{}

Real
HydrostaticBC::computeQpResidual()
{
  Real traction[3];
  for(int i=0;i<3;i++){
    traction[i]=-_pressure*_normals[_qp](i);
  }
  return -_test[_i][_qp]*traction[_component];
}
