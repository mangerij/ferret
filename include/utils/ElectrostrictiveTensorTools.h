#ifndef ELECTROSTRICTIVETENSORTOOLS_H
#define ELECTROSTRICTIVETENSORTOOLS_H

class RankFourTensor;

namespace ElectrostrictiveTensorTools
{

  Real computeProduct(const RankFourTensor & Cijkl, const RankFourTensor & Qmnkl);

  /// Sum over (j,l) q_ijkl*v(j)*w(l)
  Real electrostrictiveProduct(unsigned int i,const RealVectorValue & v, unsigned int k, const RealVectorValue & w);

  /// Sum over l q_ijkl*v(j)
  Real electrostrictiveProduct(unsigned int i, const RealVectorValue & v,unsigned int k, const unsigned int l);


}

#endif //ELECTROSTRICTIVETENSORTOOLS_H
