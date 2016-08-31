/**
 * @file   RotatedBulkEnergyDerivativeSixth.C
 * @author J. Mangeri <john.mangeri@uconn.edu
 *
 */

#include "RotatedBulkEnergyDerivativeSixth.h"
#include<cmath>

template<>
InputParameters validParams<RotatedBulkEnergyDerivativeSixth>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotatedBulkEnergyDerivativeSixth::RotatedBulkEnergyDerivativeSixth(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _Euler_angles(getParam<Real>("euler_angle_1"),
                 getParam<Real>("euler_angle_2"),
                 getParam<Real>("euler_angle_3")),
   _alpha1(getParam<Real>("alpha1")),
   _alpha11(getParam<Real>("alpha11")),
   _alpha12(getParam<Real>("alpha12")),
   _alpha111(getParam<Real>("alpha111")),
   _alpha112(getParam<Real>("alpha112")),
   _alpha123(getParam<Real>("alpha123")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotatedBulkEnergyDerivativeSixth::computeQpResidual()
{
  const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1) ? _polar_y: _polar_z;
  const VariableValue & _polar_j = (_component == 0) ? _polar_y : (_component == 1) ? _polar_z: _polar_x;
  const VariableValue & _polar_k = (_component == 0) ? _polar_z : (_component == 1) ? _polar_x: _polar_y;
  if(_component == 0)
    {
      RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
      RotationTensor R(_Euler_angles);
      RealVectorValue Rp = (R(0,0) * p(0) + R(1,0)* p(1) + R(2,0) * p(2), R(0,1) * p(0) + R(1,1)* p(1) + R(2,1) * p(2), R(0,2) * p(0) + R(1,2)* p(1) + R(2,2) * p(2));
      Real Rbulk = 0.0;
      Rbulk += ((2.0 * _alpha1 * Rp(0) + 4.0 * _alpha11 * std::pow(Rp(0), 3.0) + 2.0 * _alpha12 * Rp(0) * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) +
	  6.0 * _alpha111 * std::pow(Rp(0), 5.0) + 4.0 * _alpha112 * std::pow(Rp(0), 3.0) * (Rp(1) * Rp(1) + Rp(2) * Rp(2)) +
	  2.0 * _alpha112 * Rp(0) * (std::pow(Rp(1), 4.0) + std::pow(Rp(2), 4.0)) + 2.0 * _alpha123 * Rp(0) * std::pow(Rp(1), 2.0) * std::pow(Rp(2), 2.0)) * _test[_i][_qp]) * std::pow(_len_scale, 3.0);
      ///  Moose::out << "\n R_bulk-"; std::cout << _component << " = " << Rbulk;
      return Rbulk;
    }
  else if(_component == 1)
    {// p(0) -> p(1), p(1) -> p(2), p(2) -> p(0)
      RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
      RotationTensor R(_Euler_angles);
      RealVectorValue Rp = (R(0,0) * p(1) + R(1,0)* p(2) + R(2,0) * p(0), R(0,1) * p(1) + R(1,1)* p(2) + R(2,1) * p(0), R(0,2) * p(1) + R(1,2)* p(2) + R(2,2) * p(0));
      Real Rbulk = 0.0;
      Rbulk += ((2.0 * _alpha1 * Rp(0) + 4.0 * _alpha11 * std::pow(Rp(0), 3.0) + 2.0 * _alpha12 * Rp(0) * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) +
	  6.0 * _alpha111 * std::pow(Rp(0), 5.0) + 4.0 * _alpha112 * std::pow(Rp(0), 3.0) * (Rp(1) * Rp(1) + Rp(2) * Rp(2)) +
	  2.0 * _alpha112 * Rp(0) * (std::pow(Rp(1), 4.0) + std::pow(Rp(2), 4.0)) + 2.0 * _alpha123 * Rp(0) * std::pow(Rp(1), 2.0) * std::pow(Rp(2), 2.0)) * _test[_i][_qp]) * std::pow(_len_scale, 3.0);
      ///  Moose::out << "\n R_bulk-"; std::cout << _component << " = " << Rbulk;
      return Rbulk;
    }
  else
    {// p(0) -> p(2), p(1) -> p(0), p(2) -> p(1)
      RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
      RotationTensor R(_Euler_angles);
      RealVectorValue Rp = (R(0,0) * p(2) + R(1,0)* p(0) + R(2,0) * p(1), R(0,1) * p(2) + R(1,1)* p(0) + R(2,1) * p(1), R(0,2) * p(2) + R(1,2)* p(0) + R(2,2) * p(1));
      Real Rbulk = 0.0;
      Rbulk += ((2.0 * _alpha1 * Rp(0) + 4.0 * _alpha11 * std::pow(Rp(0), 3.0) + 2.0 * _alpha12 * Rp(0) * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) +
	  6.0 * _alpha111 * std::pow(Rp(0), 5.0) + 4.0 * _alpha112 * std::pow(Rp(0), 3.0) * (Rp(1) * Rp(1) + Rp(2) * Rp(2)) +
	  2.0 * _alpha112 * Rp(0) * (std::pow(Rp(1), 4.0) + std::pow(Rp(2), 4.0)) + 2.0 * _alpha123 * Rp(0) * std::pow(Rp(1), 2.0) * std::pow(Rp(2), 2.0)) * _test[_i][_qp]) * std::pow(_len_scale, 3.0);
      ///  Moose::out << "\n R_bulk-"; std::cout << _component << " = " << Rbulk;
      return Rbulk;
    }
}

Real
RotatedBulkEnergyDerivativeSixth::computeQpJacobian()
{
  if(_component == 0)
    {
      const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1) ? _polar_y: _polar_z;
      const VariableValue & _polar_j = (_component == 0) ? _polar_y : (_component == 1) ? _polar_z: _polar_x;
      const VariableValue & _polar_k = (_component == 0) ? _polar_z : (_component == 1) ? _polar_x: _polar_y;
      RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
      RotationTensor R(_Euler_angles);
      RealVectorValue Rp = (R(0,0) * p(0) + R(1,0)* p(1) + R(2,0) * p(2), R(0,1) * p(0) + R(1,1)* p(1) + R(2,1) * p(2), R(0,2) * p(0) + R(1,2)* p(1) + R(2,2) * p(2));
      Real Rbulk = 0.0;
      return (2.0 * _alpha1 + 12.0 * _alpha11 * std::pow(Rp(0), 2) +
	  2.0 * _alpha12 * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) + 30.0 * _alpha111 * std::pow(Rp(0), 4.0) +
	  12.0 * _alpha112 * std::pow(Rp(0), 2.0) * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) + 2.0 * _alpha112 * (std::pow(Rp(1), 4.0) + std::pow(Rp(2), 4.0)) +
	  2.0 * _alpha123 * std::pow(Rp(1), 2.0) * std::pow(Rp(2), 2.0)
  ) * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
    }
  else if(_component == 1)
    {// p(0) -> p(1), p(1) -> p(2), p(2) -> p(0)
      const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1) ? _polar_y: _polar_z;
      const VariableValue & _polar_j = (_component == 0) ? _polar_y : (_component == 1) ? _polar_z: _polar_x;
      const VariableValue & _polar_k = (_component == 0) ? _polar_z : (_component == 1) ? _polar_x: _polar_y;

      RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);

      RotationTensor R(_Euler_angles);

      RealVectorValue Rp = (R(0,0) * p(1) + R(1,0)* p(2) + R(2,0) * p(0), R(0,1) * p(1) + R(1,1)* p(2) + R(2,1) * p(0), R(0,2) * p(1) + R(1,2)* p(2) + R(2,2) * p(0));
      return (2.0 * _alpha1 + 12.0 * _alpha11 * std::pow(Rp(0), 2) +
	  2.0 * _alpha12 * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) + 30.0 * _alpha111 * std::pow(Rp(0), 4.0) +
	  12.0 * _alpha112 * std::pow(Rp(0), 2.0) * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) + 2.0 * _alpha112 * (std::pow(Rp(1), 4.0) + std::pow(Rp(2), 4.0)) +
	  2.0 * _alpha123 * std::pow(Rp(1), 2.0) * std::pow(Rp(2), 2.0)  ) * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
    }
  else
    {// p(0) -> p(2), p(1) -> p(0), p(2) -> p(1)
      const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1) ? _polar_y: _polar_z;
      const VariableValue & _polar_j = (_component == 0) ? _polar_y : (_component == 1) ? _polar_z: _polar_x;
      const VariableValue & _polar_k = (_component == 0) ? _polar_z : (_component == 1) ? _polar_x: _polar_y;
      RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
      RotationTensor R(_Euler_angles);
      RealVectorValue Rp = (R(0,0) * p(2) + R(1,0)* p(0) + R(2,0) * p(1), R(0,1) * p(2) + R(1,1)* p(0) + R(2,1) * p(1), R(0,2) * p(2) + R(1,2)* p(0) + R(2,2) * p(1));
      return (2.0 * _alpha1 + 12.0 * _alpha11 * std::pow(Rp(0), 2) +
	  2.0 * _alpha12 * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) + 30.0 * _alpha111 * std::pow(Rp(0), 4.0) +
	  12.0 * _alpha112 * std::pow(Rp(0), 2.0) * (std::pow(Rp(1), 2.0) + std::pow(Rp(2), 2.0)) + 2.0 * _alpha112 * (std::pow(Rp(1), 4.0) + std::pow(Rp(2), 4.0)) +
	  2.0 * _alpha123 * std::pow(Rp(1), 2.0) * std::pow(Rp(2), 2.0)
  ) * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
    }
}

Real
RotatedBulkEnergyDerivativeSixth::computeQpOffDiagJacobian(unsigned int jvar)
{

  mooseAssert(jvar != variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var)
    {
    if(_component == 0)
      {
        Real r;
        const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1)? _polar_y: _polar_z;
        const VariableValue & _polar_j = (jvar == _polar_x_var)? _polar_x : (jvar == _polar_y_var)? _polar_y: _polar_z;
        const VariableValue & _polar_k = ((_component == 0 && jvar == _polar_y_var) || (_component == 1 && jvar == _polar_x_var) )? _polar_z : ( (_component == 0 && jvar == _polar_z_var) || (_component == 2 && jvar == _polar_x_var))? _polar_y: _polar_x;
        RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
        RotationTensor R(_Euler_angles);
        RealVectorValue Rp = (R(0,0) * p(0) + R(1,0)* p(1) + R(2,0) * p(2), R(0,1) * p(0) + R(1,1)* p(1) + R(2,1) * p(2), R(0,2) * p(0) + R(1,2)* p(1) + R(2,2) * p(2));
        r = (4.0 * _alpha12 * Rp(0) * Rp(1) + 8.0 * _alpha112 * std::pow(Rp(0), 3.0) * Rp(1)
        + 8.0 *_alpha112 * Rp(0) * std::pow(Rp(1), 3.0) + 4.0 * _alpha123 * Rp(0) * Rp(1) * std::pow(Rp(2), 2.0));
        return r * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
      }
    else if(_component == 1)
      {// p(0) -> p(1), p(1) -> p(2), p(2) -> p(0)
        Real r;
        const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1)? _polar_y: _polar_z;
        const VariableValue & _polar_j = (jvar == _polar_x_var)? _polar_x : (jvar == _polar_y_var)? _polar_y: _polar_z;
        const VariableValue & _polar_k = ((_component == 0 && jvar == _polar_y_var) || (_component == 1 && jvar == _polar_x_var) )? _polar_z : ( (_component == 0 && jvar == _polar_z_var) || (_component == 2 && jvar == _polar_x_var))? _polar_y: _polar_x;
        RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
        RotationTensor R(_Euler_angles);
        RealVectorValue Rp = (R(0,0) * p(1) + R(1,0)* p(2) + R(2,0) * p(0), R(0,1) * p(1) + R(1,1)* p(2) + R(2,1) * p(0), R(0,2) * p(1) + R(1,2)* p(2) + R(2,2) * p(0));
        r = (4.0 * _alpha12 * Rp(0) * Rp(1) + 8.0 * _alpha112 * std::pow(Rp(0), 3.0) * Rp(1)
        + 8.0 *_alpha112 * Rp(0) * std::pow(Rp(1), 3.0) + 4.0 * _alpha123 * Rp(0) * Rp(1) * std::pow(Rp(2), 2.0));
        return r * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
      }
    else
      {// p(0) -> p(2), p(1) -> p(0), p(2) -> p(1)
        Real r;
        const VariableValue & _polar_i = (_component == 0) ? _polar_x : (_component == 1)? _polar_y: _polar_z;
        const VariableValue & _polar_j = (jvar == _polar_x_var)? _polar_x : (jvar == _polar_y_var)? _polar_y: _polar_z;
        const VariableValue & _polar_k = ((_component == 0 && jvar == _polar_y_var) || (_component == 1 && jvar == _polar_x_var) )? _polar_z : ( (_component == 0 && jvar == _polar_z_var) || (_component == 2 && jvar == _polar_x_var))? _polar_y: _polar_x;
        RealVectorValue p = (_polar_i[_qp], _polar_j[_qp], _polar_k[_qp]);
        RotationTensor R(_Euler_angles);
        RealVectorValue Rp = (R(0,0) * p(2) + R(1,0)* p(0) + R(2,0) * p(1), R(0,1) * p(2) + R(1,1)* p(0) + R(2,1) * p(1), R(0,2) * p(2) + R(1,2)* p(0) + R(2,2) * p(1));
        r = (4.0 * _alpha12 * Rp(0) * Rp(1) + 8.0 * _alpha112 * std::pow(Rp(0), 3.0) * Rp(1)
        + 8.0 *_alpha112 * Rp(0) * std::pow(Rp(1), 3.0) + 4.0 * _alpha123 * Rp(0) * Rp(1) * std::pow(Rp(2), 2.0));
        return r * _test[_i][_qp] * _phi[_j][_qp] * std::pow(_len_scale, 3.0);
      }
    }
  else
    return 0.0;
}
