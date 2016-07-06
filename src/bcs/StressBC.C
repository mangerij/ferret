/*************************************************************************
*  Stolen back from Hedgehog
****************************************************************/

#include "StressBC.h"

template<>
InputParameters validParams<StressBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredParam<int>("component","Which component(0 for x, 1 for y, 2 for z) in traction is used");

    params.addParam<std::vector<Real> >("boundary_stress", "Boundary stress: s11, s22, s33, s23, s13, s12");

    params.addCoupledVar("boundary_stress_vars", "Variable names for the: s11, s22, s33");

    params.addParam<bool>("convert_to_gpa", false, "Convert the input amounts to GPa from Pa");

    return params;
}

StressBC::StressBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _stress_vector(getParam<std::vector<Real> >("boundary_stress")),
    _Jacobian_mult(getMaterialProperty<RankFourTensor>("Jacobian_mult")),
    _component(getParam<int>("component")),
    _convert_to_gpa(getParam<bool>("convert_to_gpa")),
    _multiplier(1)
{
  if(isParamValid("boundary_stress"))
    _boundary_stress.fillFromInputVector(_stress_vector);
  else if(isCoupled("boundary_stress_vars"))
  {
    int n = coupledComponents("boundary_stress_vars");

    if(n > 3)
      mooseError("Can only take the diagonals as coupled values right now!");

    _boundary_stress_vars.resize(n);

    for (unsigned int i=0; i<_boundary_stress_vars.size(); ++i)
      _boundary_stress_vars[i] = &coupledValue("boundary_stress_vars", i);
  }
  else
    mooseError("Must provide something to StressBC!");

  if(_convert_to_gpa)
    _multiplier = 1e-9;
}

Real
StressBC::computeQpResidual()
{
  // If nothing was coupled this will be a no-op
  for(unsigned int i=0; i<_boundary_stress_vars.size(); i++)
    _boundary_stress(i, i) = _multiplier*(*_boundary_stress_vars[i])[_qp];

  return -_test[_i][_qp]*(_boundary_stress.row(_component)*_normals[_qp]);
}

Real
StressBC::computeQpJacobian()
{
  return 0;
}
