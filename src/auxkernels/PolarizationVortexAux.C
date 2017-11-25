/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "PolarizationVortexAux.h"
/// This kernel defines the polarization of a planar vortex in the xy plane characterized by its real core size c and complex (a pair of reals) location a = (a_x, a_y) in the plane.
/// The vortex evolves artificially by setting one of the three parameters a_x, a_y or c to _t -- the current time in a transient executioner.
/// Each instance of the kernel operates on one of the polarization components P_x, P_y, or P_z, renamed P_0, but couples to the other two, renamed P_1 and P_2, and 
/// ordered so that P_0, P_1, P_2 form a right-handed triple.

template<>
InputParameters validParams<PolarizationVortexAux>()
{
  InputParameters params = validParams<AuxKernel>(); params += validParams<FerretBase>();
  params.addParam<std::string>("class","PolarizationVortexAux","Class name");
  params.addRequiredParam<unsigned int>("i", "An integer corresponding to the direction the variable this kernel acts on (0 for x, 1 for y, 2 for z).");
  params.addParam<Real>("R", 1, "Magnetic dot radius");
  params.addParam<Real>("L", 1, "Magnetic dot thickness");
  params.addParam<Real>("a_x",0,"Vortex location along x");
  params.addParam<Real>("a_y", 0, "Vortex location along y");
  params.addParam<Real>("c", 1, "Vortex size");
  params.addParam<std::string>("p", "a_x", "Vortex evolution parameter: a_x, a_y or c.");
  ///params.addRequiredCoupledVar("Pnormal", "The components of polarization in the plane complementary to the i-th  direction");
  return params;
}

PolarizationVortexAux::PolarizationVortexAux(const InputParameters & parameters)
  : FerretBase(parameters), AuxKernel(parameters),
   _i(getParam<unsigned int>("i")),
   _p(getParam<std::string>("p")),
   _a_x(getParam<Real>("a_x")),
   _a_y(getParam<Real>("a_y")),
   _c(getParam<Real>("c")),
   _R(getParam<Real>("R")),
  _L(getParam<Real>("L"))
    ///    _P_1(coupledValue("Pnormal",0)),
    ///    _P_2(coupledValue("Pnormal",1))
{
  if(debug("PolarizationVortexAux")) {
      libMesh::out << "PolarizationVortexAux:\n";
      libMesh::out << "_i = " << _i << ", _a_x = " << _a_x << ", _a_y = " << _a_y << ", _c = " << _c << ", _p = " << _p << "\n";
  }
  {
    std::stringstream err;
    err << "Expected positive R, received " << _R;
    mooseAssert(_R >  0, err.str());
  }
  {
    std::stringstream err;
    err << "Expected positive L, received " << _L;
    mooseAssert(_L >  0, err.str());
  }
}


Real
PolarizationVortexAux::computeValue()
{
  Real  Pval, P = _u[_qp];
  Point xyzpoint = _current_node[_qp]; /// Assume this is a nodal kernel; recll that _current_node is a pointer
  Real  x = xyzpoint(0)/_R, y = xyzpoint(1)/_R;
  Real  a_x, a_y, c;
   if(_p == "a_x") {
    a_x = _t;
    a_y = _a_y;
    c = _c;
  }
  else if(_p == "a_y") {
    a_x = _a_x;
    a_y = _t;
    c   = _c;
  }
  else {
    a_x = _a_x;
    a_y = _a_y;
    c = _t;
  }
   {
     std::stringstream err;
     err << "Expected nonnegative c, received " << c;
     mooseAssert(c >= 0, err.str());
   }
   if(debug("computeValue")){
      libMesh::out << "PolarizationVortexAux::computeValue\n";
     libMesh::out << "a_x = " << a_x << ", a_y = " << a_y << ", c = " << c << ", i = " << _i << ", x = " << x << ", y = " << y << ", P = " << P << "\n";
   }
   Real qx,qy, q,qq; /// qx+iqy is an auxiliary complex quantity = 1 - a^*z, where z = x + iy, a = a_x + i*a_y; qq is its absolute value squared, q is the absolute value
   qx = 1.0 - (a_x*x+a_y*y);
   qy = y*a_x - x*a_y;
   qq = qx*qx+qy*qy;
   q = sqrt(qq);
   Real xa,ya,Ra, Rc; /// shift relative to a, distance from the center and core radius
   xa = x - a_x;
   ya = y - a_y;
   Ra = sqrt(xa*xa + ya*ya);
   Rc = c*sqrt(qq);
   if (debug("computeValue")) {
     libMesh::out << "qq = " << qq << ", q = " << q << ", Ra = " << Ra << ", Rc = " << Rc << "\n";
   }
   
  Real  wx, wy, ww;
  if(Ra < Rc) {
    wx = (xa*qx+ya*qy)/(q*Rc);
    wy = (-xa*qy+ya*qx)/(q*Rc);
  }
  else {
    wx = (xa*qx+ya*qy)/(q*Ra);
    wy = (-xa*qy+ya*qx)/(q*Ra);
  }
  ww = wx*wx+wy*wy;
  if (debug("computeValue")) {
    libMesh::out << "ww = " << ww << ", wx = " << wx << ", wy = " << wy << "\n";
  }
  switch(_i) { /// x, y or z?
  case 0:
    Pval = 2.0*wx/(1.0+ww);
    break;
  case 1:
    Pval = 2.0*wy/(1.0+ww);
    break;
  case 2:
    Pval = (1.0-ww)/(1+ww);
    break;
  }
  if (debug("computeValue")) {
    libMesh::out << "Pval = " << Pval << "\n";
  }
  return Pval;
}

