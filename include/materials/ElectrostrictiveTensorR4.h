/**
 * @file   ElectrostrictiveTensorR4.h
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Mon Nov 25 11:30:40 2013
 * @brief  ElectrostrictiveTensorR4 is a rank four tensor; the entries are indexed by i,j,k,l equal to 0,1,2
 *         It holds 81 entries, however, only 36 entries are independent due to the constrait symmetry constraint q_ijkl=q_jikl=q_ijlk
 */
#ifndef ELECTROSTRICTIVETENSORR4_H
#define ELECTROSTRICTIVETENSORR4_H
#include "RankFourTensor.h"
#include "ElasticityTensorR4.h"
#include "ElectrostrictiveCoefficientR4.h"
class ElectrostrictiveTensorR4 : public RankFourTensor

{
public:

  void computeProduct(const ElasticityTensorR4& Cijkl, const ElasticityTensorR4& Qmnkl);
  Real electrostrictiveProduct(unsigned int i,const RealVectorValue& v, unsigned int k, const RealVectorValue& w)const; //Sum over (j,l) q_ijkl*v(j)*w(l)
  Real electrostrictiveProduct(unsigned int i, const RealVectorValue& v,unsigned int k, const unsigned int l)const; //Sum over l q_ijkl*v(j)
};

#endif //ELECTROSTRICTIVETENSORR4_H
