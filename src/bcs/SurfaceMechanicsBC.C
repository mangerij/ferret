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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "SurfaceMechanicsBC.h"
#include "IntegratedBC.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "MooseMesh.h"

registerMooseObject("FerretApp", SurfaceMechanicsBC);

InputParameters SurfaceMechanicsBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<std::vector<Real> >("Cs_ijkl", "Surface elastic tensor,  C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222"); //read (at most) 9 components of the surface elastic tensor
  params.addRequiredParam<RealVectorValue>("S_k", "Surface euler angle vector,  S_1, S_2, S_3"); //read (at most) 9 components of the surface elastic tensor
  params.addRequiredParam<Real>("taus", "Intrinsic surface stress"); // Intrinsic surface tension is a 2x2 rank-2 tensor
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

SurfaceMechanicsBC::SurfaceMechanicsBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _component(getParam<unsigned int>("component")),
  _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain")),
  _Csijkl_vector(getParam<std::vector<Real> >("Cs_ijkl")),
  _S_k_vector(getParam<RealVectorValue>("S_k")), //fill method for RotationTensor only accepts RealVectorValue
  _taus(getParam<Real>("taus"))
{
}

RankTwoTensor
SurfaceMechanicsBC::computeQpProjection() //Compute the projection operator at a given quadrature point:
{
  RankTwoTensor projection;

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = i; j < 3; ++j)
    {
      if (i == j)
        projection(i, i) = 1.0 - _normals[_qp](i) * _normals[_qp](i);
      else
      {
       	projection(j, i) = -_normals[_qp](i) * _normals[_qp](j);
       	projection(i, j) = -_normals[_qp](i) * _normals[_qp](j);
      }
    }
  return projection;
}

RankTwoTensor
SurfaceMechanicsBC::computeQpRankTwoRotation(unsigned int tensor_component) /// compute rank-2 and rank-4 tensors used to write surface strain/stress and elastic constants in local representation
{
  RealVectorValue T1;
  RealVectorValue T2;
  RankTwoTensor a;

  T1(0) = _normals[_qp](2) - _normals[_qp](1); // Form a vector perpendicular to the surface normal by taking cross-product of (1,1,1) with _normals[_qp]
  T1(1) = _normals[_qp](0) - _normals[_qp](2);
  T1(2) = _normals[_qp](1) - _normals[_qp](0);
  if (T1 * T1 < 1.e-12) // need to check that the norm of T1 is not zero
  {
    T1(0) = -1.0; //_normal[_qp] is parallel to (1,1,1) so choose a tangent vector orthogonal to (1,1,1)
    T1(1) = 0.0;
    T1(2) = 1.0;
  }
  T1 = T1 / std::sqrt(T1 * T1); //normalize

  T2(0) = _normals[_qp](1) * T1(2) - _normals[_qp](2) * T1(1); // if A and B are normalized, then A \times B = C will also be normalized
  T2(1) = _normals[_qp](2) * T1(0) - _normals[_qp](0) * T1(2);
  T2(2) = _normals[_qp](0) * T1(1) - _normals[_qp](1) * T1(0);

  for (unsigned int i = 0; i < 3; ++i) //compute rank-two and rank-four tensors
      {
      for (unsigned int j = 0; j < 3; ++j)
        {
          if (tensor_component == 11)
            a(i,j) = T1(i) * T1(j);
          else if (tensor_component == 22)
            a(i,j) = T2(i) * T2(j);
        }
      }
  return a;
}

RankFourTensor
SurfaceMechanicsBC::computeQpRankFourRotation(unsigned int tensor_component) /// compute rank-2 and rank-4 tensors used to write surface strain/stress and elastic constants in local representation
{
  RankFourTensor a;
  RealVectorValue T1;
  RealVectorValue T2;

  T1(0) = _normals[_qp](2) - _normals[_qp](1); // Form a vector perpendicular to the surface normal by taking cross-product of (1,1,1) with _normals[_qp]
  T1(1) = _normals[_qp](0) - _normals[_qp](2);
  T1(2) = _normals[_qp](1) - _normals[_qp](0);
  if (T1 * T1 < 1.e-12) // need to check that the norm of T1 is not zero
  {
    T1(0) = -1.0; //_normal[_qp] is parallel to (1,1,1) so choose a tangent vector orthogonal to (1,1,1)
    T1(1) = 0.0;
    T1(2) = 1.0;
  }
  T1 = T1 / std::sqrt(T1 * T1); //normalize

  T2(0) = _normals[_qp](1) * T1(2) - _normals[_qp](2) * T1(1); // if A and B are normalized, then A \times B = C will also be normalized
  T2(1) = _normals[_qp](2) * T1(0) - _normals[_qp](0) * T1(2);
  T2(2) = _normals[_qp](0) * T1(1) - _normals[_qp](1) * T1(0);

  for (unsigned int i = 0; i < 3; ++i) //compute rank-two and rank-four tensors
      {
      for (unsigned int j = 0; j < 3; ++j)
        {
        for (unsigned int k = 0; k < 3; ++k)
          {
            for (unsigned int l = 0; l < 3; ++l)
              {
                if (tensor_component == 11)
                  a(i, j, k, l) = T1(i) * T1(j) * T1(k) * T1(l);
                else if (tensor_component == 22)
                  a(i, j, k, l) = T2(i) * T2(j) * T2(k) * T2(l);
                else if (tensor_component == 12)
                  a(i, j, k, l) = T1(i) * T1(j) * T2(k) * T2(l) + T2(i) * T2(j) * T1(k) * T1(l);
                else if (tensor_component == 33)
                  a(i, j, k, l) = T1(i) * T2(j) * T1(k) * T2(l) + T2(i) * T1(j) * T2(k) * T1(l) + T1(i) * T2(j) * T2(k) * T1(l) + T2(i) * T1(j) * T1(k) * T2(l);
              }
           }
         }
      }
  return a;
}

Real
SurfaceMechanicsBC::computeQpResidual()
{
  RankTwoTensor projection;

  RankTwoTensor tp11;
  RankTwoTensor tp22;

  RankFourTensor t11;
  RankFourTensor t22;
  RankFourTensor t12;
  RankFourTensor t33;

  projection = computeQpProjection(); //compute projection tensor

  tp11 = computeQpRankTwoRotation(11); //compute tpij tensors
  tp22 = computeQpRankTwoRotation(22);

  t11 = computeQpRankFourRotation(11); //compute tij tensors
  t22 = computeQpRankFourRotation(22);
  t12 = computeQpRankFourRotation(12);
  t33 = computeQpRankFourRotation(33);

  RankFourTensor Csijkl;
  Csijkl.surfaceFillFromInputVector(_Csijkl_vector);  //fill surface stiffness tensor
  RotationTensor R(_S_k_vector); //construct rotation matrix
  Csijkl.rotate(R); //Euler rotate the surface stiffness tensor in global representation

  RankTwoTensor surface_strain = projection * (_elastic_strain[_qp] * projection); // this is the surface strain in global representation
  RankFourTensor surface_elastic_constants = t11 * (Csijkl(0, 0, 0, 0) - _taus) + t22 * (Csijkl(1, 1, 1, 1) - _taus) + t12 * (Csijkl(0, 0, 1, 1) - _taus) + t33 * (Csijkl(0, 1, 0, 1) - _taus); //surface elastic constants in global representation

  RankTwoTensor surface_stress = surface_elastic_constants * surface_strain * projection; // temporary tensor surface elasticity * surface_strain in global representation
  RankTwoTensor surface_tau = (tp11 + tp22) * projection; //temporary tensor for intrinsic surface tension

  Real Rs = 0.0; //surface mechanics residual contribution
  for (unsigned int j = 0; j < 3; ++j)
    for (unsigned int k = 0; k < 3; ++k)
	    Rs += _taus * surface_tau(j, k) * _grad_test[_i][_qp](k) * projection(j, _component) + projection(j, _component) * _grad_test[_i][_qp](k) * surface_stress(j, k); //note the second term is zero for the initial nonlinear and linear step.
  return Rs;
}
