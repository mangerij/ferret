/****************************************************************/
/* Stress BC:                                                   */
/*     This BC is intended only for testing purpose.            */
/*     Mathematically, Stress BC is not well-posed.             */
/*     Anyway, it reads six functions:                          */
/*     stress_xx, stress_xy, stress_yy,                         */
/*     stress_yz, stress_zx, stress_zz                          */
/*     Then, it computes the traction by multiplying            */
/*     stress tensor with the normal to get the traction.       */
/*     The traction is used for the problem.                    */
/*                                                              */
/****************************************************************/


#include "TractionBC.h" 
template<>
InputParameters validParams<TractionBC>()
{
    InputParameters params = validParams<IntegratedBC>();
  // Here we are adding a parameter that will be extracted from the input file by the Parser
  params.addRequiredParam<int>("component","Which component(0 for x,1 for y, 2 for z) in traction is used");
  params.addRequiredParam<Real>("stress_xx", "stress_xx"); 
  params.addRequiredParam<Real>("stress_xy", "stress_xy"); 
  params.addRequiredParam<Real>("stress_yy", "stress_yy"); 
  params.addRequiredParam<Real>("stress_yz", "stress_yz"); 
  params.addRequiredParam<Real>("stress_zx", "stress_zx");
  params.addRequiredParam<Real>("stress_zz", "stress_zz"); 
  //params.addRequiredParam<Real>("traction", "traction")
  //params.addRequiredParam<RealVectorValue>("vector_value",RealVectorValue(), "traction this BC should act on");
  return params;
}

TractionBC::TractionBC(const std::string & name, InputParameters parameters) :
  IntegratedBC(name, parameters),
  _component(getParam<int>("component")),
  _stress_xx(getParam<Real>("stress_xx")),
  _stress_xy(getParam<Real>("stress_xy")),
  _stress_yy(getParam<Real>("stress_yy")),
  _stress_yz(getParam<Real>("stress_yz")),
  _stress_zx(getParam<Real>("stress_zx")),
  _stress_zz(getParam<Real>("stress_zz"))
  //_traction(getParam<Real>("traction"))
{}

Real
TractionBC::computeQpResidual()
{
  Real _values[3][3];
  Real _traction[3];
  _values[0][0]=_stress_xx; _values[0][1]=_stress_xy;_values[0][2]=_stress_zx;
  _values[1][0]=_stress_xy; _values[1][1]=_stress_yy;_values[1][2]=_stress_yz;
  _values[2][0]=_stress_zx; _values[2][1]=_stress_yz;_values[2][2]=_stress_zz;
  for(int i=0;i<3;i++){
   _traction[i]=0;
    for(int j=0;j<3;j++)
     _traction[i]=_traction[i]+_values[i][j]*_normals[_qp](j); 
  }
  return -_test[_i][_qp]*_traction[_component]; //be careful with the sign
  //  return 0.0;
}
