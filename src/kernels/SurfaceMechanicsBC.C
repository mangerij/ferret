/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "SurfaceMechanicsBC.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

#include "RotationTensor.h"

#include "MooseMesh.h"


template<>
InputParameters validParams<SurfaceMechanicsBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  // First read the (at most) nine components of the surface elastic tensor
  // C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
  // These will have to be put into a 3x3x3x3 rank-4 tensor with appropriate zero padding
  params.addRequiredParam<std::vector<Real> >("Cs_ijkl", "Surface elastic tensor,  C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222");
  // The intrinsic surface tension is a 2x2 rank-2 tensor that has to be expanded to a 3x3
  params.addRequiredParam<Real>("taus", "Intrinsic surface stress");
  params.addParam<Real>("surface_euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("surface_euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("surface_euler_angle_3", 0.0, "Euler angle in direction 3");
  params.addRequiredCoupledVar("disp_x", "The x displacement");
  params.addCoupledVar("disp_y", "The y displacement");
  params.addCoupledVar("disp_z", "The z displacement");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  return params;

}

SurfaceMechanicsBC::SurfaceMechanicsBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _dim(_mesh.dimension()),
  _component(getParam<unsigned int>("component")),
  _Csijkl_vector(getParam<std::vector<Real> >("Cs_ijkl")),
  _taus(getParam<Real>("taus")),
  _surface_Euler_angle_1(getParam<Real>("surface_euler_angle_1")),
  _surface_Euler_angle_2(getParam<Real>("surface_euler_angle_2")),
  _surface_Euler_angle_3(getParam<Real>("surface_euler_angle_3")),
  _surface_Euler_angles(_surface_Euler_angle_1, _surface_Euler_angle_2, _surface_Euler_angle_3),
  _grad_disp_x(coupledGradient("disp_x")),
  _grad_disp_y(_dim >= 2 ? coupledGradient("disp_y") : _grad_zero),
  _grad_disp_z(_dim == 3 ? coupledGradient("disp_z") : _grad_zero)
{
}

void
SurfaceMechanicsBC::computeQpProjection()
{
  RankTwoTensor projection;

  //Compute the projection operator at a given quadrature point:
  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=i; j<3; ++j)
    {
      if (i == j)
        projection(i, i) = 1.0 - _normals[_qp](i) * _normals[_qp](i);
      else
      {
       	projection(j, i) = -_normals[_qp](i) * _normals[_qp](j);
       	projection(i, j) = -_normals[_qp](i) * _normals[_qp](j);
      }
    }
}

void
SurfaceMechanicsBC::computeQpRotation()
/// compute rank-4 tensors used to write surface elastic constants in local representation
/// First form a vector perpendicular to the surface normal by taking cross-product of (1,1,1) with _normals[_qp]
{
  RealVectorValue tangent_1;
  RealVectorValue tangent_2;

  RankFourTensor t11;
  RankFourTensor t22;
  RankFourTensor t12;
  RankFourTensor t33;

  RankTwoTensor tp11;
  RankTwoTensor tp22;

  tangent_1(0) = _normals[_qp](2) - _normals[_qp](1);
  tangent_1(1) = _normals[_qp](0) - _normals[_qp](2);
  tangent_1(2) = _normals[_qp](1) - _normals[_qp](0);
  if (tangent_1 * tangent_1 < 1.e-12) // need to check that the norm of _tangent_1 is not zero
    {
    tangent_1(0) = -1.0;     //_normal[_qp] is parallel to (1,1,1) so choose a tangent vector orthogonal to (1,1,1)
    tangent_1(1) = 0.0;
    tangent_1(2) = 1.0;
    }
    tangent_1 = tangent_1/sqrt(tangent_1 * tangent_1);
    tangent_2(0) = _normals[_qp](1) * tangent_1(2) - _normals[_qp](2) * tangent_1(1);
    tangent_2(1) = _normals[_qp](2) * tangent_1(0) - _normals[_qp](0) * tangent_1(2);
    tangent_2(2) = _normals[_qp](0) * tangent_1(1) - _normals[_qp](1) * tangent_1(0);

    for (unsigned int i = 0; i < 3; ++i)  //compute rank-two and rank-four tensors
      {
      for (unsigned int j = 0; j < 3; ++j)
        {
         tp11(i,j) = tangent_1(i) * tangent_1(j);
         tp22(i,j) = tangent_2(i) * tangent_2(j);
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
  RankTwoTensor grad_tensor_1(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);
  RankTwoTensor temp2;
  RankTwoTensor temp4;

  RankTwoTensor surface_strain;
  RankTwoTensor surface_stress;
  RankTwoTensor projection;

  RankTwoTensor tp11;
  RankTwoTensor tp22;

  RankFourTensor temp_surface;
  RankFourTensor t11;
  RankFourTensor t22;
  RankFourTensor t12;
  RankFourTensor t33;

  RankFourTensor Csijkl;
  Csijkl.surfaceFillFromInputVector (_Csijkl_vector);
  RotationTensor R(_surface_Euler_angles);
  Csijkl.rotate(R);

  Real C0000 = Csijkl(1, 1, 1, 1) - _taus; //apply the tension
  Real C1111 = Csijkl(2, 2, 2, 2) - _taus;
  Real C0011 = Csijkl(1, 1, 2, 2) + _taus;
  Real C0101 = Csijkl(1, 2, 1, 2) - _taus;

  temp2 = (grad_tensor_1 + grad_tensor_1.transpose())/2.0;
  // Project onto surface
  computeQpProjection();
  computeQpRotation();

  surface_strain = projection * (temp2 * projection); // this is the surface strain in global representation
  temp_surface = t11 * C0000 + t22 * C1111 + t12 * C0011 + t33 * C0101; //surface elastic constants in global representation

  surface_stress = temp_surface * surface_strain; //this is now surface stress (projected onto tangent plane) using global representation
  RankTwoTensor grad_tensor_2(_grad_test[_i][_qp], _grad_test[_i][_qp], _grad_test[_i][_qp]);
  RankTwoTensor temp_tau;
  // temp4 = _projection*(grad_tensor_2*_projection); // projected tensor of gradient of test functions
  temp4 = surface_stress * projection; // temporary tensor surface elasticity
  temp_tau = (tp11 + tp22) * projection; // temporary tensor for intrinsic surface stress
  Real temp = 0.0;
  unsigned int cp;
  unsigned int ip;
  unsigned int jp;
  cp = _component;
  for (unsigned int i = 0; i < 3; ++i)
    {
    ip=i;
      for (unsigned int j = 0; j < 3; ++j)
        {
	      jp=j;
	      // temp+=_surface_tau(_component,i)*temp4(_component,i)+_projection(i,_component)*_grad_test[_i][_qp](j)*temp4(i,j);
	      // temp+=_surface_tau.(cp,ip)*temp2.(cp,ip)+_projection.(ip,cp)*_grad_test[_i][_qp](j)*temp4.(ip,jp);
	      // temp+=_taus*((_tp11.(cp,ip)+_tp22.(cp,ip))*temp2.(cp,ip))+_projection.(ip,cp)*_grad_test[_i][_qp](j)*temp4.(ip,jp);
	      temp += _taus * temp_tau(ip, jp) * _grad_test[_i][_qp](j) * projection(ip, cp) + projection(ip,cp) * _grad_test[_i][_qp](j) * temp4(ip,jp);
        }
    }
    return temp;
}
