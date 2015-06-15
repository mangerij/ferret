/**
 * @file   ElectrostrictiveTensorR4.C
 * @author S. Gu <sgu@anl.gov>
 * @modifed J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15 2015
 * @brief  ElectrostrictiveTensorR4 is a rank four tensor; the entries are indexed by i,j,k,l equal to 0,1,2
 *         It holds 81 entries, however, only 36 entries are independent due to the symmetry constraint q_ijkl=q_jikl=q_ijlk
 */

#include "ElectrostrictiveTensorR4.h"

void
ElectrostrictiveTensorR4::computeProduct(const ElasticityTensorR4& Cijkl, const ElasticityTensorR4& Qmnkl)
{
  for(unsigned int i=0;i<3;++i)
    for(unsigned int j=0;j<3;++j)
      for(unsigned int k=0;k<3;++k)
        for(unsigned int l=0;l<3;++l)
        {
	         Real sum=0.0;
	         for(unsigned int m=0;m<3;++m)
	           for(unsigned int n=0;n<3;++n)
           {
	             sum += Cijkl(i,j,m,n)*Qmnkl(m,n,k,l);
	         }

        _vals[i][j][k][l] = 2*sum;
    }
}


Real ElectrostrictiveTensorR4::electrostrictiveProduct(unsigned int i, const RealVectorValue& v,unsigned int k, const unsigned int l)const
{
  //Sum over j q_ijkl*v(j)
  Real sum=0.0;
  for(unsigned int j=0;j<3;++j)
  {
    sum+=(*this)(i,j,k,l)*v(j);
  }
  return sum;
}
Real ElectrostrictiveTensorR4::electrostrictiveProduct(unsigned int i,const RealVectorValue& v,unsigned int k, const RealVectorValue& w)const
{
  //Sum over (j,l) q_ijkl*v(j)*w(l)
  Real sum=0.0;
  for(unsigned int j=0;j<3; ++j)
    for(unsigned int l=0;l<3;++l){
      sum+=(*this)(i,j,k,l)*v(j)*w(l);
    }
  return sum;
}
