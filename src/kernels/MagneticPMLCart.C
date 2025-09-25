/*
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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "MagneticPMLCart.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticPMLCart);

InputParameters MagneticPMLCart::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution to Laplacian in a stretched region");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. 0 for x, 1 for y, 2 for z).");
  /*   params.addParam<Real>("deltasxminus",1.e03, "scaled thickness of pml region (dimensionful)");
   params.addParam<Real>("deltapxminus",1., "pole distance");
   params.addParam<Real>("deltawxminus",1.,"distance to be scaled (dimensionfull");
   params.addParam<Real>("x0pmlminus",-1., "position at which scaling will start");*/
   params.addParam<Real>("deltasyminus", 1.e03, "scaled thickness of pml region");
   params.addParam<Real>("deltapyminus",1., "pole distance");
   params.addParam<Real>("deltawyminus",1.,"distance to be scaled (dimensionfull");
   params.addParam<Real>("y0pmlminus",-1., "position at which scaling will start");
   /*   params.addParam<Real>("deltaszminus",1.e03, "sclaed thickness of pml region");
   params.addParam<Real>("deltapzminus",1., "pole distance");
   params.addParam<Real>("deltawzminus",1.,"distance to be scaled (dimensionfull");
   params.addParam<Real>("z0pmlminus",-1., "position at which scaling will start");
   params.addParam<Real>("deltasxplus",1.e03, "scaled thickness of pml region");
   params.addParam<Real>("deltapxplus",1., "pole distance");
   params.addParam<Real>("deltawxplus",1.,"distance to be scaled (dimensionfull");
   params.addParam<Real>("x0pmlplus",1., "position at which scaling will start");
   params.addParam<Real>("deltasyplus", 1.e03, "distance to be scaled along y");
   params.addParam<Real>("deltapyplus",1., "pole distance");
   params.addParam<Real>("deltawyplus",1.,"distance to be scaled (dimensionfull");
   params.addParam<Real>("y0pmlplus",1., "position at which scaling will start");
   params.addParam<Real>("deltaszplus",1.e03, "scaled thickness of pml region");
   params.addParam<Real>("deltapzplus",1., "pole distance");
   params.addParam<Real>("deltawzplus",1.,"distance to be scaled (dimensionfull");
   params.addParam<Real>("z0pmlplus",1., "position at which scaling will start");*/
  return params;
}

MagneticPMLCart::MagneticPMLCart(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
    _permittivity(getMaterialProperty<Real>("permittivity")),
  /*    _deltasxminus(getParam<Real>("deltasxminus")),
    _deltapxminus(getParam<Real>("deltapxminus")),
    _deltawxminus(getParam<Real>("deltawxminus")),
    _x0pmlminus(getParam<Real>("x0pmlminus")),*/
    _deltasyminus(getParam<Real>("deltasyminus")),
    _deltapyminus(getParam<Real>("deltapyminus")),
    _deltawyminus(getParam<Real>("deltawyminus")),
    _y0pmlminus(getParam<Real>("y0pmlminus"))
    /*    _deltaszminus(getParam<Real>("deltaszminus")),
    _deltapzminus(getParam<Real>("deltapzminus")),
    _deltawzminus(getParam<Real>("deltawxminus")),
    _z0pmlminus(getParam<Real>("z0pmlminus")),
    _deltasxplus(getParam<Real>("deltasxplus")),
    _deltapxplus(getParam<Real>("deltapxplus")),
    _deltawxplus(getParam<Real>("deltawxminus")),
    _x0pmlplus(getParam<Real>("x0pmlplus")),
    _deltasyplus(getParam<Real>("deltasyplus")),
    _deltapyplus(getParam<Real>("deltapyplus")),
    _deltawyplus(getParam<Real>("deltawxminus")),
    _y0pmlplus(getParam<Real>("y0pmlplus")),
    _deltaszplus(getParam<Real>("deltaszplus")),
    _deltapzplus(getParam<Real>("deltapzplus")),
    _deltawzplus(getParam<Real>("deltawxminus")),
    _z0pmlplus(getParam<Real>("z0pmlplus"))*/
{
}

Real
MagneticPMLCart::computeQpResidual()
{
  /*
  if (_component == 0)
    {
      if (_q_point[_qp](0) > _x0pmlplus)
	  {
	  const Real gamma = (_deltasxplus + _deltapxplus)/_deltasxplus;
          const Real xi = (_q_point[_qp](0)-_x0pmlplus)/_deltawxplus;
          const Real dudx =(1.+1./(gamma-xi))*_deltapxplus/(_deltawxplus*(gamma-xi));
          return _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp]/dudx;
	  }
	  else if (_q_point[_qp](0) < _x0pmlminus)
	  {
	  const Real gamma = (_deltasxminus + _deltapxminus)/_deltasxminus;
          const Real xi = -(_q_point[_qp](0)-_x0pmlminus)/_deltawxminus;
          const Real dudx = (1.+1./(gamma-xi))*_deltapxminus/(_deltawxminus*(gamma-xi));
          return _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp]/dudx;
	  }
          else
	    return 0.0;

    }
  else
  */
  if (_component == 1)
  {
    /*
    if (_q_point[_qp](1) > _y0pmlplus)
         {
	 const Real gamma = (_deltasyplus + _deltapyplus)/_deltasyplus;
         const Real xi = (_q_point[_qp](1)-_y0pmlplus)/_deltawyplus;
         const Real dudx = 1.+_deltapyplus/_deltawyplus*(xi/(gamma-xi))*((2.*gamma-1.)/(gamma-1.));
         return _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp]/dudx;
         }
	 else
    */
	 if  (_q_point[_qp](1) < _y0pmlminus)
	 {
	   const Real gamma = (_deltasyminus + _deltapyminus)/_deltasyminus;
	   const Real xi = -(_q_point[_qp](1)-_y0pmlminus)/_deltawyminus;
//	   const Real dudx = 1.+_deltapyminus/_deltawyminus*(xi/(gamma-xi))*((2.*gamma-xi)/(gamma-xi));
           const Real dudx = 1.+_deltapyminus/_deltawyminus*(1./(gamma-xi)*(1.+xi/(gamma-xi))-1./gamma);
	   return _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp]/dudx;
          }
	   else
	   return 0.0;
    }

  /*else
   if (_component == 2)
  {
   if (_q_point[_qp](2) > _z0pmlplus)
         {
	 const Real gamma = (_deltaszplus + _deltapzplus)/_deltaszplus;
         const Real xi = (_q_point[_qp](2)-_z0pmlplus)/_deltawzplus;
         const Real dudx = (1.+1./(gamma-xi))*_deltapzplus/(_deltawzplus*(gamma-xi));
         return _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp]/dudx;
         }
	 else if  (_q_point[_qp](2) < _z0pmlminus)
	 {
	 const Real gamma = (_deltaszminus + _deltapzminus)/_deltaszminus;
         const Real xi = -(_q_point[_qp](2)-_z0pmlminus)/_deltawzminus;
         const Real dudx = (1.+1./(gamma-xi))*_deltapzminus/(_deltawzminus*(gamma-xi));
         return _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp]/dudx;
          }
         else
	   return 0.0;
  }
  */
   else
     return 0.0;
}

Real
MagneticPMLCart::computeQpJacobian()
{
  /*    if (_component == 0)
    {
      if (_q_point[_qp](0) > _x0pmlplus)
	  {
	  const Real gamma = (_deltasxplus + _deltapxplus)/_deltasxplus;
          const Real xi = (_q_point[_qp](0)-_x0pmlplus)/_deltawxplus;
          const Real dudx = (1.+1./(gamma-xi))*_deltapxplus/(_deltawxplus*(gamma-xi));
          return _permittivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp]/dudx;
	  }
	  else if (_q_point[_qp](0) < _x0pmlminus)
	  {
	  const Real gamma = (_deltasxminus + _deltapxminus)/_deltasxminus;
          const Real xi = -(_q_point[_qp](0)-_x0pmlminus)/_deltawxminus;
          const Real dudx =  (1.+1./(gamma-xi))*_deltapxplus/(_deltawxplus*(gamma-xi));
          return _permittivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp]/dudx;
	  }
          else
	    return 0.0;

    }
  else
*/
  if (_component == 1)
  {
    /*  if (_q_point[_qp](1) > _y0pmlplus)
         {
	 const Real gamma = (_deltasyplus + _deltapyplus)/_deltasyplus;
         const Real xi = (_q_point[_qp](1)-_y0pmlplus)/_deltawyplus;
         const Real dudx =  1.+_deltapyplus/_deltawyplus*(xi/(gamma-xi))*((2.*gamma-1.)/(gamma-1.));
         return _permittivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp]/dudx;
         }
	 else */
	 if  (_q_point[_qp](1) < _y0pmlminus)
	 {
      	   const Real gamma = (_deltasyminus + _deltapyminus)/_deltasyminus;
           const Real xi = -(_q_point[_qp](1)-_y0pmlminus)/_deltawyminus;
//           const Real  dudx = 1.+_deltapyminus/_deltawyminus*(xi/(gamma-xi))*((2.*gamma-xi)/(gamma-xi));
           const Real dudx = 1.+_deltapyminus/_deltawyminus*(1./(gamma-xi)*(1.+xi/(gamma-xi))-1./gamma);
	   return _permittivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp]/dudx;
          }
         else
	   return 0.0;
  }
/*
  else if (_component == 2)
  {
   if (_q_point[_qp](2) > _z0pmlplus)
         {
	 const Real gamma = (_deltaszplus + _deltapzplus)/_deltaszplus;
         const Real xi = (_q_point[_qp](2)-_z0pmlplus)/_deltawzplus;
         const Real dudx = (1.+1./(gamma-xi))*_deltapzplus/(_deltawzplus*(gamma-xi));
         return _permittivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp]/dudx;
         }
	 else if  (_q_point[_qp](2) < _z0pmlminus)
	 {
      	 const Real gamma = (_deltaszminus + _deltapzminus)/_deltaszminus;
         const Real xi = -(_q_point[_qp](2)-_z0pmlminus)/_deltawzminus;
         const Real dudx = (1.+1./(gamma-xi))*_deltapzminus/(_deltawzminus*(gamma-xi));
	 return _permittivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp]/dudx;
          }
         else
	   return 0.0;
	   }*/

  else
    return 0.0;
    }

