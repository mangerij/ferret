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


#include "StressFunctionBC.h"
#include "Function.h" 
template<>
InputParameters validParams<StressFunctionBC>()
{
    InputParameters params = validParams<IntegratedBC>();
  // Here we are adding a parameter that will be extracted from the input file by the Parser
  params.addRequiredParam<int>("component","Which component(0 for x,1 for y, 2 for z) in traction is used");
  params.addRequiredParam<FunctionName>("stress_xx", "stress_xx(function)"); 
  params.addRequiredParam<FunctionName>("stress_xy", "stress_xy(function)"); 
  params.addRequiredParam<FunctionName>("stress_yy", "stress_yy(function)"); 
  params.addRequiredParam<FunctionName>("stress_yz", "stress_yz(function)"); 
  params.addRequiredParam<FunctionName>("stress_zx", "stress_zx(function)");
  params.addRequiredParam<FunctionName>("stress_zz", "stress_zz(function)"); 
  return params;
}

StressFunctionBC::StressFunctionBC(const std::string & name, InputParameters parameters) :
  IntegratedBC(name, parameters),
  _component(getParam<int>("component")),
  _stress_xx(getFunction("stress_xx")),
  _stress_xy(getFunction("stress_xy")),
  _stress_yy(getFunction("stress_yy")),
  _stress_yz(getFunction("stress_yz")),
  _stress_zx(getFunction("stress_zx")),
  _stress_zz(getFunction("stress_zz"))
{
}

Real
StressFunctionBC::computeQpResidual()
{
  Real _values[3][3];
  Real _traction[3];
  _values[0][0]=_stress_xx.value(_t,_q_point[_qp]); _values[0][1]=_stress_xy.value(_t,_q_point[_qp]);_values[0][2]=_stress_zx.value(_t,_q_point[_qp]);
  _values[1][0]=_stress_xy.value(_t,_q_point[_qp]); _values[1][1]=_stress_yy.value(_t,_q_point[_qp]);_values[1][2]=_stress_yz.value(_t,_q_point[_qp]);
  _values[2][0]=_stress_zx.value(_t,_q_point[_qp]); _values[2][1]=_stress_yz.value(_t,_q_point[_qp]);_values[2][2]=_stress_zz.value(_t,_q_point[_qp]);
  for(int i=0;i<3;i++){
    _traction[i]=0;
    for(int j=0;j<3;j++)
      _traction[i]=_traction[i]+_values[i][j]*_normals[_qp](j); 
  }
  return -_test[_i][_qp]*_traction[_component]; //be careful with the sign
}
