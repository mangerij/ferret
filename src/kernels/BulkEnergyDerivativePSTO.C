
/**
 * @file   BulkEnergyDerivativePSTO.C
 * @author J. Mangeri <john.mangeri@uconn.edu> and S. Churchill <steven.churchill@uconn.edu>
 *
 */

#include "BulkEnergyDerivativePSTO.h"
#include<cmath>

template<>
InputParameters validParams<BulkEnergyDerivativePSTO>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha3", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha4", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha5", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x3", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x4", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x5", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x6", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("epsilon", "Constant value for strain");
  params.addRequiredParam<Real>("T", "Temperature");
  params.addRequiredParam<Real>("Tc", "Critical Temperature");
  return params;
}

BulkEnergyDerivativePSTO::BulkEnergyDerivativePSTO(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _alpha1(getParam<Real>("alpha1")),
   _alpha2(getParam<Real>("alpha2")),
   _alpha3(getParam<Real>("alpha3")),
   _alpha4(getParam<Real>("alpha4")),
   _alpha5(getParam<Real>("alpha5")),
   _x1(getParam<Real>("x1")),
   _x2(getParam<Real>("x2")),
   _x3(getParam<Real>("x3")),
   _x4(getParam<Real>("x4")),
   _x5(getParam<Real>("x5")),
   _x6(getParam<Real>("x6")),
   _epsilon(getParam<Real>("epsilon")),
   _T(getParam<Real>("T")),
   _Tc(getParam<Real>("Tc"))
{
  std::cout<<"_alpha1 ="<<_alpha1<<"\n";
  std::cout<<"_alpha2 ="<<_alpha2<<"\n";
  std::cout<<"_alpha3 ="<<_alpha3<<"\n";
  std::cout<<"_alpha4 ="<<_alpha4<<"\n";
  std::cout<<"_alpha5 ="<<_alpha5<<"\n";
  std::cout<<"_x1 ="<<_x1<<"\n";
  std::cout<<"_x2 ="<<_x2<<"\n";
  std::cout<<"_x3 ="<<_x3<<"\n";
  std::cout<<"_x4 ="<<_x4<<"\n";
  std::cout<<"_x5 ="<<_x5<<"\n";
  std::cout<<"_x6 ="<<_x6<<"\n";
  std::cout<<"_epsilon ="<<_epsilon<<"\n";
  std::cout<<"_T ="<<_T<<"\n";
  std::cout<<"_Tc ="<<_Tc<<"\n";

}

Real
BulkEnergyDerivativePSTO::computeQpResidual()
{
  if (_component == 0)
  {
   Real Efree = 0.0;
   Efree += (_alpha1 * (_T - _Tc) * (2.0 * _polar_x[_qp]) + _alpha2 * 4.0 * std::pow(_polar_x[_qp], 3.0) + 
_alpha3 * 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0) + 6.0 * _alpha4 * std::pow(_polar_x[_qp], 5.0) + 
_alpha5 * (4.0 * std::pow(_polar_x[_qp], 3.0) * std::pow(_polar_y[_qp], 2.0) + 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 4.0)) + 
(_x1 * ( 2.0 * _polar_x[_qp]) + _x2 * (4.0 * std::pow(_polar_x[_qp], 3.0)) + 
_x3 * 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0)) * _epsilon + (_x4 * (2.0 * _polar_x[_qp]) + 
_x5 * ( 4.0 * std::pow(_polar_x[_qp], 3.0)) + _x6 * 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0)) * std::pow(_epsilon, 2.0))  * _test[_i][_qp];    
    return Efree;
  }
  else if (_component == 1)
  {
   Real Efree = 0.0;
   Efree += (_alpha1 * (_T - _Tc) * (2.0 * _polar_y[_qp] ) + _alpha2 * 4.0 * std::pow(_polar_y[_qp], 3.0) + 
_alpha3 * 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0) + 6.0 * _alpha4 *std::pow(_polar_y[_qp], 5.0) + 
_alpha5 * (4.0 * std::pow(_polar_y[_qp], 3.0) * std::pow(_polar_x[_qp], 2.0) + 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 4.0)) + 
(_x1 * ( 2.0 * _polar_y[_qp]) + _x2 * (4.0 * std::pow(_polar_y[_qp], 3.0)) + 
_x3 * 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0)) * _epsilon + (_x4 * (2.0 * _polar_y[_qp]) + 
_x5 * ( 4.0 * std::pow(_polar_y[_qp], 3.0)) + _x6 * 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0)) * std::pow(_epsilon, 2.0))  * _test[_i][_qp];
   return Efree;
  }
  else 
  {
    return 0;
  }
}

Real
BulkEnergyDerivativePSTO::computeQpJacobian()
{
  if (_component == 0)
  {
    Real EfreeJac = 0.0;
    EfreeJac += (_alpha1 * (_T - _Tc) * 2.0 + _alpha2 * 12.0 * std::pow(_polar_x[_qp], 2.0) + 
_alpha3 * 2.0 * std::pow(_polar_y[_qp], 2.0) + 30.0 * _alpha4 * std::pow(_polar_x[_qp], 4.0) + 
_alpha5 * (12.0 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + 2.0 * std::pow(_polar_y[_qp], 4.0)) + 
(_x1 * 2.0 + _x2 * 12.0 * std::pow(_polar_x[_qp], 2.0) + 
_x3 * 2.0 * std::pow(_polar_y[_qp], 2.0)) * _epsilon + (_x4 * 2.0 + 
_x5 * 12.0 * std::pow(_polar_x[_qp], 2.0) + _x6 * 2.0 * std::pow(_polar_y[_qp], 2.0)) * std::pow(_epsilon, 2.0)) * _phi[_j][_qp] * _test[_i][_qp];    
    return EfreeJac;
  }
  else if (_component == 1)
  {
    Real EfreeJac = 0.0;
    EfreeJac += (_alpha1 * (_T - _Tc) * 2.0 + _alpha2 * 12.0 * std::pow(_polar_y[_qp], 2.0) + 
_alpha3 * 2.0 * std::pow(_polar_x[_qp], 2.0) + 30.0 * _alpha4 * std::pow(_polar_y[_qp], 4.0) + 
_alpha5 * (12.0 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_x[_qp], 2.0) + 2.0 * std::pow(_polar_x[_qp], 4.0)) + 
(_x1 * 2.0 + _x2 * 12.0 * std::pow(_polar_y[_qp], 2.0) + 
_x3 * 2.0 * std::pow(_polar_x[_qp], 2.0)) * _epsilon + (_x4 * 2.0 + 
_x5 * 12.0 * std::pow(_polar_y[_qp], 2.0) + _x6 * 2.0 * std::pow(_polar_x[_qp], 2.0)) * std::pow(_epsilon, 2.0)) * _phi[_j][_qp] * _test[_i][_qp];
    return EfreeJac;
  }
  else 
  {
    return 0;
  }
}

Real
BulkEnergyDerivativePSTO::computeQpOffDiagJacobian()
{
  if (_component == 0 or 1)
  {
    Real EfreeJac = 0.0;
    EfreeJac += (_alpha3 * 4.0 * _polar_x[_qp] * _polar_y[_qp] + 
_alpha5 * (8.0 * std::pow(_polar_x[_qp], 3.0) * _polar_y[_qp] + 8.0 * _polar_x[_qp] * std::pow(+polar_y[_qp], 3.0)) + 
(_x3 * 4.0 * _polar_x[_qp] * _polar_y[_qp]) * _epsilon +
(_x6 * 4.0 * _polar_x[_qp] * _polar_y[_qp]) * std::pow(_epsilon, 2.0)) * _test[_i][_qp] * _phi[_j][_qp]; 
    return EfreeJac;
  }

  else 
  {
    return 0;
  }
}




