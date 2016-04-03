/****************************************************************/
/* SurfaceMechanicsBC                                          */
/*     Includes surface stress contributions to the free energy.*/
/*                                                              */
/*    IntegrateBC: the independent surface elasticity constants */
/*    and Euler angles for each surface. A residual is formed   */
/*    as an integrated boundary condition, as IntegratedBC      */
/*    has access to surface integrals and surface normals       */
/*                                                              */
/****************************************************************/
/* Note that the function surfacefillFromVector was added       */
/* to RankTwoTensor and RankFourTensor in tensor_mechanics      */
/* O. Heinonen and A. Jokisaari                                 */

#include "Kernel.h"
#include "SurfaceMechanicsBC.h"
#include "IntegratedBC.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

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
    params.addRequiredCoupledVar("disp_y", "The y displacement");
    params.addCoupledVar("disp_z", "The z displacement");
    params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
    return params;

    //    InputParameters params2 = validParams<Kernel>();
    //    params2.addRequiredParam<unsigned int>("component", "an integer corresponding to the direction the variable in the kernel acts on");
    //    return params2;
}

SurfaceMechanicsBC::SurfaceMechanicsBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _dim(_mesh.dimension()),
  _component(getParam<unsigned int>("component")),
  _surface_euler_angle_1(getParam<Real>("surface_euler_angle_1")),
  _surface_euler_angle_2(getParam<Real>("surface_euler_angle_2")),
  _surface_euler_angle_3(getParam<Real>("surface_euler_angle_3")),
  _Csijkl_vector(getParam<std::vector<Real> >("Cs_ijkl")),
  _taus(getParam<Real>("taus")),
  _Csijkl(),
  _surface_euler_angles(_surface_euler_angle_1, _surface_euler_angle_2, _surface_euler_angle_3),
  _grad_disp_x(coupledGradient("disp_x")),
  _grad_disp_y(coupledGradient("disp_y")),
  _grad_disp_z(_dim == 3 ? coupledGradient("disp_z") : _grad_zero)
{
  _Csijkl.surfaceFillFromInputVector (_Csijkl_vector);
  RotationTensor R(_surface_euler_angles);
  _Csijkl.rotate(R);
}

void
SurfaceMechanicsBC::computeQpProjection()
{
  //Compute the projection operator at a given quadrature point:
  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=i; j<3; ++j)
    {
      if (i == j)
	    {
        _projection(i, i) = 1 - _normals[_qp](i) * _normals[_qp](i);
      }
      else
      {
       	_projection(j, i) = -_normals[_qp](i) * _normals[_qp](j);
       	_projection(i, j) = -_normals[_qp](i) * _normals[_qp](j);
      }
    }
}

void
SurfaceMechanicsBC :: computeQpRotation()
// compute rank-4 tensors used to write surface elastic constants in local representation
// First form a vector perpendicular to the surface normal by taking cross-product of (1,1,1) with _normals[_qp]
{
  _tangent_1(0) = _normals[_qp](2) - _normals[_qp](1);
  _tangent_1(1) = _normals[_qp](0) - _normals[_qp](2);
  _tangent_1(2) = _normals[_qp](1) - _normals[_qp](0);
// need to check that the norm of _tangent_1 is not zero
  if (_tangent_1 * _tangent_1 < 1.e-12)
    {
    //_normal[_qp] is parallel to (1,1,1) so choose a tangent vector orthogonal to (1,1,1)
    _tangent_1(0) = -1.;
    _tangent_1(1) = 0.;
    _tangent_1(2) = 1.;
    }
    _tangent_1 = _tangent_1/sqrt(_tangent_1 * _tangent_1);
    _tangent_2(0) = _normals[_qp](1) * _tangent_1(2) - _normals[_qp](2) * _tangent_1(1);
    _tangent_2(1) = _normals[_qp](2) * _tangent_1(0) - _normals[_qp](0) * _tangent_1(2);
    _tangent_2(2) = _normals[_qp](0) * _tangent_1(1) - _normals[_qp](1) * _tangent_1(0);
//compute rank-two and rank-four tensors
    for (unsigned int i = 0; i < 3; ++i)
      {
      for (unsigned int j = 0; j < 3; ++j)
        {
         _tp11(i,j) = _tangent_1(i) * _tangent_1(j);
         _tp22(i,j) = _tangent_2(i) * _tangent_2(j);
        for (unsigned int k = 0; k < 3; ++k)
	        {
          for (unsigned int l = 0; l < 3; ++l)
	          {
	          _t11(i, j, k, l) = _tangent_1(i) * _tangent_1(j) * _tangent_1(k) * _tangent_1(l);
	          _t22(i, j, k, l) = _tangent_2(i) * _tangent_2(j) * _tangent_2(k) * _tangent_2(l);
	          _t12(i, j, k, l) = _tangent_1(i) * _tangent_1(j) * _tangent_2(k) * _tangent_2(l) + _tangent_2(i) * _tangent_2(j) * _tangent_1(k) * _tangent_1(l);
	          _t33(i, j, k, l) = _tangent_1(i) * _tangent_2(j) * _tangent_1(k) * _tangent_2(l) + _tangent_2(i) * _tangent_1(j) * _tangent_2(k) * _tangent_1(l) + _tangent_1(i) * _tangent_2(j) * _tangent_2(k) * _tangent_1(l) + _tangent_2(i) * _tangent_1(j) * _tangent_1(k) * _tangent_2(l);
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
  RankFourTensor temp_surface;
  C0000 = _Csijkl(1, 1, 1, 1) - _taus;
  C1111 = _Csijkl(2, 2, 2, 2) - _taus;
  C0011 = _Csijkl(1, 1, 2, 2) + _taus;
  C0101 = _Csijkl(1, 2, 1, 2) - _taus;
  //_surface_strain=(grad_tensor_1+grad_tensor_1.transpose())/2.0;
  temp2 = (grad_tensor_1 + grad_tensor_1.transpose())/2.0;
  //Project onto surface
  computeQpProjection();
  computeQpRotation();
  //  temp2=_projection*(_surface_strain*_projection); // this is the surface strain in global representation
  _surface_strain = _projection * (temp2 * _projection); // this is the surface strain in global representation
  temp_surface = _t11 * C0000 + _t22 * C1111 + _t12 * C0011 + _t33*C0101; //surface elastic constants in global representation
  //_surface_stress=temp_surface*temp2; //this is now surface stress (projected onto tangent plane) using global representation
  _surface_stress = temp_surface * _surface_strain; //this is now surface stress (projected onto tangent plane) using global representation
  RankTwoTensor grad_tensor_2(_grad_test[_i][_qp], _grad_test[_i][_qp], _grad_test[_i][_qp]);
  RankTwoTensor temp_tau;
  //temp4=_projection*(grad_tensor_2*_projection);// projected tensor of gradient of test functions
   temp4 = _surface_stress * _projection;// temporary tensor surface elasticity
   temp_tau = (_tp11 + _tp22) * _projection;//temporary tensor for intrinsic surface stress
   Real temp;
   temp = 0.0;
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
	  //temp+=_surface_tau(_component,i)*temp4(_component,i)+_projection(i,_component)*_grad_test[_i][_qp](j)*temp4(i,j);
	  //temp+=_surface_tau.(cp,ip)*temp2.(cp,ip)+_projection.(ip,cp)*_grad_test[_i][_qp](j)*temp4.(ip,jp);
	  //temp+=_taus*((_tp11.(cp,ip)+_tp22.(cp,ip))*temp2.(cp,ip))+_projection.(ip,cp)*_grad_test[_i][_qp](j)*temp4.(ip,jp);
	  temp += _taus * temp_tau(ip, jp) * _grad_test[_i][_qp](j) * _projection(ip, cp) + _projection(ip,cp) * _grad_test[_i][_qp](j) * temp4(ip,jp);
        }
    }
   	return temp;
     //Why no computeQpOffDiagJacobian() or computeQpJacobian() ?
}
