#include "PolarizationSurfaceCharge.h"

template<>
InputParameters validParams<PolarizationSurfaceCharge>()
{
  InputParameters params = validParams<IntegratedBC>(); params += validParams<FerretBase>();
  params.addParam<std::string>("class", "PolarizationSurfacCharge", "Class name");
  params.addRequiredCoupledVar("P", "Polarization");
  params.addParam<bool>("J_polarization",false, "Whether to compute the Jacobian block with respect to polarization");
  return params;
}

PolarizationSurfaceCharge::PolarizationSurfaceCharge(const InputParameters & parameters) :
  FerretBase(parameters), IntegratedBC(parameters),
  _P_x(coupledValue("P", 0)),
  _P_y(coupledValue("P", 1)),
  _P_z(coupledValue("P", 2)),
  _J_polarization(getParam<bool>("J_polarization"))
{
}

Real
PolarizationSurfaceCharge::computeQpResidual()
{
  RealVectorValue P(_P_x[_qp], _P_y[_qp], _P_z[_qp]);
  RealVectorValue n(_normals[_qp]);

  Real charge = n*P;
  Real res = _test[_i][_qp]*charge;
  if (debug("computeQpResidual"))
  {
    libMesh::out << "PolarizationSurfaceCharge::computeQpResidual:\n";
    libMesh::out << "_qp = " << _qp << ", _qpoint[_qp] = " << _q_point[_qp] << ", n = " << n << ", P = " << P << ", charge = " << charge << ", res = " << res  << "\n";
  }
  return res;
}
