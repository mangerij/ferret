/****************************************************************/
/* HydrostaticDirichelt BC                                      */
/****************************************************************/


#include "HydrostaticDirichletBC.h"
template<>
InputParameters validParams<HydrostaticDirichletBC>()
{
    InputParameters params = validParams<NodalBC>();
  // Here we are adding a parameter that will be extracted from the input file by the Parser
  params.addRequiredParam<Real>("center_x", "x coordinate of center of sphere");
  params.addRequiredParam<Real>("center_y", "y coordinate of center of sphere");
  params.addRequiredParam<Real>("center_z", "z coordinate of center of sphere");
  params.addRequiredParam<Real>("displacement","total inward displacement");
  params.addRequiredParam<int>("component","component of the displacement");

return params;
}

HydrostaticDirichletBC::HydrostaticDirichletBC(const InputParameters & parameters) :
  NodalBC(parameters),
  _center_x(getParam<Real>("center_x")),
  _center_y(getParam<Real>("center_y")),
  _center_z(getParam<Real>("center_z")),
  _displacement(getParam<Real>("displacement")),
  _component(getParam<int>("component"))
{}

Real
HydrostaticDirichletBC::computeQpResidual()
{
	Point _center;
	_center(0)=_center_x;
        _center(1)=_center_y;
        _center(2)=_center_z;
	Point normal_direction=(_center-*_current_node);
	Point normal=normal_direction/normal_direction.norm();

  return _u[_qp]+normal(_component)*_displacement; //be careful with the sign
  //  return 0.0;
}
