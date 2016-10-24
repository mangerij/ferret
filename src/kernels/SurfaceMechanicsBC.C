/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "SurfaceMechanicsBC.h"
#include "IntegratedBC.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<SurfaceMechanicsBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  //params.addRequiredCoupledVar("disp_x", "The x displacement");
  //params.addRequiredCoupledVar("disp_y", "The y displacement");
  //params.addCoupledVar("disp_z", "The z displacement");
  params.addRequiredParam<std::vector<Real> >("Cs_ijkl", "Surface elastic tensor,  C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222"); //read (at most) 9 components of the surface elastic tensor
  params.addRequiredParam<Real>("taus", "Intrinsic surface stress"); // Intrinsic surface tension is a 2x2 rank-2 tensor
  params.addParam<Real>("surface_euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("surface_euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("surface_euler_angle_3", 0.0, "Euler angle in direction 3");
  //params.addRequiredCoupledVar("displacements", "The string of displacements suitable for the problem statement");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

SurfaceMechanicsBC::SurfaceMechanicsBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _component(getParam<unsigned int>("component")),
  //_grad_disp_x(coupledGradient("disp_x")),
  //_grad_disp_y(coupledGradient("disp_y")),
  _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain")),
  _surface_euler_angle_1(getParam<Real>("surface_euler_angle_1")),
  _surface_euler_angle_2(getParam<Real>("surface_euler_angle_2")),
  _surface_euler_angle_3(getParam<Real>("surface_euler_angle_3")),
  _Csijkl_vector(getParam<std::vector<Real> >("Cs_ijkl")),
  _taus(getParam<Real>("taus")),
  _Csijkl(),
  _surface_euler_angles(_surface_euler_angle_1, _surface_euler_angle_2, _surface_euler_angle_3)
{
  _Csijkl.surfaceFillFromInputVector (_Csijkl_vector);
  RotationTensor R(_surface_euler_angles);
  _Csijkl.rotate(R);
}

void
SurfaceMechanicsBC::computeQpProjection()
{
  //Compute the projection operator at a given quadrature point:
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = i; j < 3; ++j)
    {
      if (i == j)
        _projection(i, i) = 1.0 - _normals[_qp](i) * _normals[_qp](i);
      else
      {
       	_projection(j, i) = -_normals[_qp](i) * _normals[_qp](j);
       	_projection(i, j) = -_normals[_qp](i) * _normals[_qp](j);
      }
    }
}

void
SurfaceMechanicsBC::computeQpRotation() /// compute rank-4 tensors used to write surface elastic constants in local representation
{
  RealVectorValue tangent_1;
  RealVectorValue tangent_2;

  RankFourTensor t11;
  RankFourTensor t22;
  RankFourTensor t12;
  RankFourTensor t33;

  RankFourTensor tp11;
  RankFourTensor tp22;

  tangent_1(0) = _normals[_qp](2) - _normals[_qp](1); // First form a vector perpendicular to the surface normal by taking cross-product of (1,1,1) with _normals[_qp]
  tangent_1(1) = _normals[_qp](0) - _normals[_qp](2);
  tangent_1(2) = _normals[_qp](1) - _normals[_qp](0);
  if (tangent_1 * tangent_1 < 1.e-12) // need to check that the norm of _tangent_1 is not zero
    {  
      tangent_1(0) = -1.0; //_normal[_qp] is parallel to (1,1,1) so choose a tangent vector orthogonal to (1,1,1)
      tangent_1(1) = 0.0;
      tangent_1(2) = 1.0;
    }
      tangent_1 = tangent_1 / sqrt(tangent_1 * tangent_1);
      tangent_2(0) = _normals[_qp](1) * tangent_1(2) - _normals[_qp](2) * tangent_1(1);
      tangent_2(1) = _normals[_qp](2) * tangent_1(0) - _normals[_qp](0) * tangent_1(2);
      tangent_2(2) = _normals[_qp](0) * tangent_1(1) - _normals[_qp](1) * tangent_1(0);

    for (unsigned int i = 0; i < 3; ++i) //compute rank-two and rank-four tensors
      {
      for (unsigned int j = 0; j < 3; ++j)
        {
         _tp11(i,j) = tangent_1(i) * tangent_1(j);
         _tp22(i,j) = tangent_2(i) * tangent_2(j);
        for (unsigned int k = 0; k < 3; ++k)
	        {
          for (unsigned int l = 0; l < 3; ++l)
	          {
	          t11(i, j, k, l) = tangent_1(i) * tangent_1(j) * tangent_1(k) * tangent_1(l);
	          t22(i, j, k, l) = tangent_2(i) * tangent_2(j) * tangent_2(k) * tangent_2(l);
	          t12(i, j, k, l) = tangent_1(i) * tangent_1(j) * tangent_2(k) * tangent_2(l) + tangent_2(i) * tangent_2(j) * tangent_1(k) * tangent_1(l);
	          t33(i, j, k, l) = tangent_1(i) * tangent_2(j) * tangent_1(k) * tangent_2(l) + tangent_2(i) * tangent_1(j) * tangent_2(k) * tangent_1(l) + tangent_1(i) * tangent_2(j) * tangent_2(k) * tangent_1(l) + tangent_2(i) * tangent_1(j) * tangent_1(k) * tangent_2(l);
            }
        }
     }
  }
}

Real
SurfaceMechanicsBC::computeQpResidual()
{
  RankFourTensor t11;
  RankFourTensor t22;
  RankFourTensor t12;
  RankFourTensor t33;

  computeQpProjection(); //compute projection tensor
  computeQpRotation(); //compute tpij and tij tensors

  RankFourTensor temp_surface;

  RankTwoTensor temp4 = _surface_stress * _projection;// temporary tensor surface elasticity

  _surface_strain = _projection * (_elastic_strain[_qp] * _projection); // this is the surface strain in global representation
  temp_surface = t11 * (_Csijkl(1, 1, 1, 1) - _taus) + t22 * (_Csijkl(2, 2, 2, 2) - _taus) + t12 * (_Csijkl(1, 1, 2, 2) + _taus) + t33 * (_Csijkl(1, 2, 1, 2) - _taus); //surface elastic constants in global representation
  _surface_stress = temp_surface * _surface_strain; //this is now surface stress (projected onto tangent plane) using global representation

  RankTwoTensor temp_tau = (_tp11 + _tp22) * _projection;//temporary tensor for intrinsic surface stress

  Real temp = 0.0;
  for (unsigned int j = 0; j < 3; ++j)
    for (unsigned int k = 0; k < 3; ++k)
	    temp += _taus * temp_tau(j, k) * _grad_test[_i][_qp](k) * _projection(j, _component) + _projection(j, _component) * _grad_test[_i][_qp](k) * temp4(j, k);
  return temp;
}
