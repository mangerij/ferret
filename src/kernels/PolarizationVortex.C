#include "PolarizationVortex.h"

template<>
InputParameters validParams<PolarizationVortex>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("P_1", "The x polarization in the oriented complementary plane");
  params.addRequiredCoupledVar("P_2", "The y polarization in the oriented complementary plane");

  return params;
}


PolarizationVortex::PolarizationVortex(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _component(getParam<unsigned int>("component")),
   _parameter(getParam<std::string>("time")),
   _a_x(getParam<Real>("a_x")),
   _a_y(getParam<Real>("a_y")),
   _c(getParam<Real>("c")),
   _P_1(coupledValue("P_1")),
   _P_2(coupledValue("P_2"))
{}

Real
PolarizationVortex::computeQpResidual()
{
  Real  Pnew, P = _u[_qp];
  Point xyzpoint = _q_point[_qp];
  Real  x = xyzpoint(0), y = xyzpoint(1);
  Real  a_x, a_y, c;
  if(_parameter == "a_x") {
    a_x = _t;
    a_y = _a_y;
    c = _c;
  }
  else if(_parameter == "a_y") {
    a_x = _a_x;
    a_y = _t;
    c   = _c;
  }
  else {
    a_x = _a_x;
    a_y = _a_y;
    c = _t;
  }
  std::stringstream err;
  err << "Expected nonnegative c, received " << c;
  mooseAssert(c >= 0, err.str());
  Real qx,qy, qq; // auxiliary complex quantity 1 - a^*z, where z = x + iy, a = a_x + i*a_y; qq is its absolute value squared
  qx = 1.0 - (a_x*x+a_y*y);
  qy = y*a_x - x*a_y;
  qq = qx*qx+qy*qy;
  Real xa,ya,R, R0; // shift relative to a, distance from the center and core radius
  xa = x - a_x;
  ya = y - a_y;
  R = sqrt(xa*xa + ya*ya);
  R0 = c*sqrt(qq);

  Real  wx, wy, ww;
  if(R < R0) {
    wx = (xa*qx+ya*qy)/(qq*c);
    wy = (xa*qy-ya*qx)/(qq*c);
  }
  else {
    wx = (xa*qx+ya*qy)/(sqrt(qq)*R);
    wy = (xa*qy-ya*qx)/(sqrt(qq)*R);
  }
  ww = wx*wx+wy*wy;

  switch(_component) { // x, y or z?
  case 0:
    Pnew = 2.0*wx/(1.0+ww);
    break;
  case 1:
    Pnew = 2.0*wy/(1.0+ww);
    break;
  case 2:
    Pnew = 0.0;
    break;
  }

  return _test[_i][_qp]*(P-Pnew);
}

Real
PolarizationVortex::computeQpJacobian()
{
  return 1.0;
}

// Real
// StressDivergence::computeQpOffDiagJacobian(unsigned int jvar)
// {
//   return 0.0;
// }
