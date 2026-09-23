#include "CubicAnisotropyCartLLG.h"

registerMooseObject("FerretApp", CubicAnisotropyCartLLG);

InputParameters
CubicAnisotropyCartLLG::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("LLG torque from cubic magnetocrystalline anisotropy, "
                             "F = Kc1 (mx^2 my^2 + my^2 mz^2 + mz^2 mx^2), single lattice.");
  params.addRequiredParam<unsigned int>("component", "Component (0, 1, 2) this kernel acts on.");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetization");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetization");
  params.addParam<Real>("g0", 1.0, "electron gyromagnetic factor");
  params.addParam<Real>("mu0", 1.0, "permeability of the vacuum");
  params.addParam<Real>("Hscale", 1.0, "scaling factor for effective fields");
  return params;
}

CubicAnisotropyCartLLG::CubicAnisotropyCartLLG(const InputParameters & parameters)
  : Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _mag_x_var(coupled("mag_x")),
    _mag_y_var(coupled("mag_y")),
    _mag_z_var(coupled("mag_z")),
    _mag_x(coupledValue("mag_x")),
    _mag_y(coupledValue("mag_y")),
    _mag_z(coupledValue("mag_z")),
    _alpha(getMaterialProperty<Real>("alpha")),
    _Kc1(getMaterialProperty<Real>("Kc1")),
    _Ms(getMaterialProperty<Real>("Ms")),
    _g0(getParam<Real>("g0")),
    _mu0(getParam<Real>("mu0")),
    _Hscale(getParam<Real>("Hscale"))
{
}

Real
CubicAnisotropyCartLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (_Kc1[_qp]*_g0*(_alpha[_qp]*_mag_x[_qp]*(pow(_mag_x[_qp], 2)*(-2*pow(_mag_y[_qp], 2) - 2*pow(_mag_z[_qp], 2)) + 2*pow(_mag_y[_qp], 4) + 2*pow(_mag_z[_qp], 4)) + _mag_y[_qp]*(-2*pow(_mag_y[_qp], 2)*_mag_z[_qp] + 2*pow(_mag_z[_qp], 3)))/(_Hscale*_Ms[_qp]*_mu0)) * _test[_i][_qp] / (1.0 + _alpha[_qp] * _alpha[_qp]);
  }
  else if (_component == 1)
  {
    return (_Kc1[_qp]*_g0*(_alpha[_qp]*(pow(_mag_x[_qp], 2)*(2*pow(_mag_x[_qp], 2)*_mag_y[_qp] - 2*pow(_mag_y[_qp], 3)) + _mag_y[_qp]*(-2*pow(_mag_y[_qp], 2)*pow(_mag_z[_qp], 2) + 2*pow(_mag_z[_qp], 4))) + _mag_x[_qp]*(2*pow(_mag_x[_qp], 2)*_mag_z[_qp] - 2*pow(_mag_z[_qp], 3)))/(_Hscale*_Ms[_qp]*_mu0)) * _test[_i][_qp] / (1.0 + _alpha[_qp] * _alpha[_qp]);
  }
  else if (_component == 2)
  {
    return (_Kc1[_qp]*_g0*(_alpha[_qp]*(pow(_mag_x[_qp], 2)*(2*pow(_mag_x[_qp], 2)*_mag_z[_qp] - 2*pow(_mag_z[_qp], 3)) + pow(_mag_y[_qp], 2)*(2*pow(_mag_y[_qp], 2)*_mag_z[_qp] - 2*pow(_mag_z[_qp], 3))) + _mag_x[_qp]*(-2*pow(_mag_x[_qp], 2)*_mag_y[_qp] + 2*pow(_mag_y[_qp], 3)))/(_Hscale*_Ms[_qp]*_mu0)) * _test[_i][_qp] / (1.0 + _alpha[_qp] * _alpha[_qp]);
  }
  else
    return 0.0;
}

Real
CubicAnisotropyCartLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_Kc1[_qp]*_alpha[_qp]*_g0*(pow(_mag_x[_qp], 2)*(-6*pow(_mag_y[_qp], 2) - 6*pow(_mag_z[_qp], 2)) + 2*pow(_mag_y[_qp], 4) + 2*pow(_mag_z[_qp], 4))/(_Hscale*_Ms[_qp]*_mu0)) * _phi[_j][_qp] * _test[_i][_qp] / (1.0 + _alpha[_qp] * _alpha[_qp]);
  }
  else if (_component == 1)
  {
    return (_Kc1[_qp]*_alpha[_qp]*_g0*(pow(_mag_x[_qp], 2)*(2*pow(_mag_x[_qp], 2) - 6*pow(_mag_y[_qp], 2)) - 6*pow(_mag_y[_qp], 2)*pow(_mag_z[_qp], 2) + 2*pow(_mag_z[_qp], 4))/(_Hscale*_Ms[_qp]*_mu0)) * _phi[_j][_qp] * _test[_i][_qp] / (1.0 + _alpha[_qp] * _alpha[_qp]);
  }
  else if (_component == 2)
  {
    return (_Kc1[_qp]*_alpha[_qp]*_g0*(pow(_mag_x[_qp], 2)*(2*pow(_mag_x[_qp], 2) - 6*pow(_mag_z[_qp], 2)) + pow(_mag_y[_qp], 2)*(2*pow(_mag_y[_qp], 2) - 6*pow(_mag_z[_qp], 2)))/(_Hscale*_Ms[_qp]*_mu0)) * _phi[_j][_qp] * _test[_i][_qp] / (1.0 + _alpha[_qp] * _alpha[_qp]);
  }
  else
    return 0.0;
}

Real
CubicAnisotropyCartLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real f = _phi[_j][_qp] * _test[_i][_qp] / (1.0 + _alpha[_qp] * _alpha[_qp]);
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
      return (_Kc1[_qp]*_g0*(_alpha[_qp]*_mag_x[_qp]*(-4*pow(_mag_x[_qp], 2)*_mag_y[_qp] + 8*pow(_mag_y[_qp], 3)) - 6*pow(_mag_y[_qp], 2)*_mag_z[_qp] + 2*pow(_mag_z[_qp], 3))/(_Hscale*_Ms[_qp]*_mu0)) * f;
    if (jvar == _mag_z_var)
      return (_Kc1[_qp]*_g0*(_alpha[_qp]*_mag_x[_qp]*(-4*pow(_mag_x[_qp], 2)*_mag_z[_qp] + 8*pow(_mag_z[_qp], 3)) + _mag_y[_qp]*(-2*pow(_mag_y[_qp], 2) + 6*pow(_mag_z[_qp], 2)))/(_Hscale*_Ms[_qp]*_mu0)) * f;
  }
  else if (_component == 1)
  {
    if (jvar == _mag_x_var)
      return (_Kc1[_qp]*_g0*(_alpha[_qp]*_mag_x[_qp]*(8*pow(_mag_x[_qp], 2)*_mag_y[_qp] - 4*pow(_mag_y[_qp], 3)) + 6*pow(_mag_x[_qp], 2)*_mag_z[_qp] - 2*pow(_mag_z[_qp], 3))/(_Hscale*_Ms[_qp]*_mu0)) * f;
    if (jvar == _mag_z_var)
      return (_Kc1[_qp]*_g0*(_alpha[_qp]*_mag_y[_qp]*(-4*pow(_mag_y[_qp], 2)*_mag_z[_qp] + 8*pow(_mag_z[_qp], 3)) + _mag_x[_qp]*(2*pow(_mag_x[_qp], 2) - 6*pow(_mag_z[_qp], 2)))/(_Hscale*_Ms[_qp]*_mu0)) * f;
  }
  else if (_component == 2)
  {
    if (jvar == _mag_x_var)
      return (_Kc1[_qp]*_g0*(_alpha[_qp]*_mag_x[_qp]*(8*pow(_mag_x[_qp], 2)*_mag_z[_qp] - 4*pow(_mag_z[_qp], 3)) - 6*pow(_mag_x[_qp], 2)*_mag_y[_qp] + 2*pow(_mag_y[_qp], 3))/(_Hscale*_Ms[_qp]*_mu0)) * f;
    if (jvar == _mag_y_var)
      return (_Kc1[_qp]*_g0*(_alpha[_qp]*_mag_y[_qp]*(8*pow(_mag_y[_qp], 2)*_mag_z[_qp] - 4*pow(_mag_z[_qp], 3)) + _mag_x[_qp]*(-2*pow(_mag_x[_qp], 2) + 6*pow(_mag_y[_qp], 2)))/(_Hscale*_Ms[_qp]*_mu0)) * f;
  }
  return 0.0;
}
