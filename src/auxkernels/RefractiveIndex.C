/**
 * @file   RefractiveIndex.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate photoelastic change to the refractive index
 * \Delta \epsilon_{ij} = - n_i^2 n_j^2 \Delta B_{ij}
 *
 * where \Delta B_{ij} = p_{ijkl} \varepsilon_{kl}
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#include "RefractiveIndex.h"
#include "RotationTensor.h"

template<>

InputParameters validParams<RefractiveIndex>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("refractive_index_bulk_ordinary", 1.0,"refractive index along the ordinary axis");
  params.addParam<Real>("refractive_index_bulk_extraordinary", 1.0,"refractive index along the extraordinary axis");
  params.addRequiredParam<unsigned int>("index_one", "A rank two tensor is being computed, need two indices");
  params.addRequiredParam<unsigned int>("index_two", "A rank two tensor is being computed, need two indices");
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}


RefractiveIndex::RefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
   _pijkl(getMaterialProperty<RankFourTensor>("photostrictive_tensor")),
   _no(getParam<Real>("refractive_index_bulk_ordinary")),
   _ne(getParam<Real>("refractive_index_bulk_extraordinary")),
   _index_one(getParam<unsigned int>("index_one")),
   _index_two(getParam<unsigned int>("index_two")),
   _Euler_angles(getParam<Real>("euler_angle_1"),
                 getParam<Real>("euler_angle_2"),
                 getParam<Real>("euler_angle_3"))
{
}

Real
RefractiveIndex::computeValue()
{
  RotationTensor R(_Euler_angles);
  // Assume that n_e is along the z-axis for now
  // note the regular birefringence is quantified by _ne - _no
  RealVectorValue n(_no, _no, _ne); 

  // Rotate the indicatrix such that it is aligned with the crystallographic direction A_i = R_{ij} A_j
  RealVectorValue nR(R(0, 0) * n(0) + R(0, 1) * n(1) + R(0, 2) * n(2), R(1, 0) * n(0) + R(1, 1) * n(1) + R(1, 2) * n(2), R(2, 0) * n(0) + R(2, 1) * n(1) + R(2, 2) * n(2));

  return -nR(_index_one)*nR(_index_one)*nR(_index_two)*nR(_index_two) * (
_pijkl[_qp](_index_one, _index_two, 0, 0) * _strain[_qp](0, 0) + 2.0 * _pijkl[_qp](_index_one, _index_two, 0, 1) * _strain[_qp](0, 1)+
2.0 * _pijkl[_qp](_index_one, _index_two, 0, 2) * _strain[_qp](0, 2) + _pijkl[_qp](_index_one, _index_two, 1, 1) * _strain[_qp](1, 1)+
2.0 * _pijkl[_qp](_index_one, _index_two, 1, 2) * _strain[_qp](1, 2) + _pijkl[_qp](_index_one, _index_two, 2, 2) * _strain[_qp](2, 2));
}


