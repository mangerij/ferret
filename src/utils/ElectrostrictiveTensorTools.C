/**
 * @file   ElectrostrictiveTensorTools.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   April 17 2016
 * @brief  ElectrostrictiveTensorR4 is a rank four tensor; the entries are indexed by i,j,k,l equal to 0,1,2
 *         It holds 81 entries, however, only 36 entries are independent due to the symmetry constraint q_ijkl=q_jikl=q_ijlk
 */

#include "MooseTypes.h"
#include "RankFourTensor.h"

namespace ElectrostrictiveTensorTools
{

RankFourTensor
computeProduct(const RankFourTensor & Cijkl, const RankFourTensor & Qmnkl)
{
  RankFourTensor result;
    // Moose::out << "\n Printing Cijkl entries:"; //keep for debugging purposes
    // for(unsigned int a = 0; a < 3; ++a)
    //   for(unsigned int b = 0; b < 3; ++b)
    //     for(unsigned int c = 0; c < 3; ++c)
    //       for(unsigned int d = 0; d < 3; ++d)
    //         {
    //           Moose::out << "\n C"; std::cout << a + 1 << b + 1 << c + 1 << d + 1; Moose::out << " = " << Cijkl(a, b, c, d) << ";";
    //         }
    // Moose::out << "\n Printing Qijkl entries:";
    // for(unsigned int r = 0; r < 3; ++r)
    //   for(unsigned int s = 0; s < 3; ++s)
    //     for(unsigned int t = 0; t < 3; ++t)
    //       for(unsigned int u = 0; u < 3; ++u)
    //         {
    //           Moose::out << "\n Q"; std::cout << r + 1 << s + 1 << t + 1 << u + 1; Moose::out << " = " << Qmnkl(r, s, t, u) << ";";
    //         }

// Moose::out << "\n Performing C_ijmn Q_mnkl contraction on all the quadrature points?";
  for(unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for(unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for(unsigned int k = 0; k < LIBMESH_DIM; ++k)
        for(unsigned int l = 0; l < LIBMESH_DIM; ++l)
        {
          Real sum = 0.0;
          for(unsigned int m = 0; m < LIBMESH_DIM; ++m)
            for(unsigned int n = 0; n < LIBMESH_DIM; ++n)
            {
              if( n != m)
              {
                //sum += 0.5 * Cijkl(i, j, m, n) * Qmnkl(m, n, k, l);
                result(i,j,k,l) += 0.5 * Cijkl(i,j,m,n) * Qmnkl(m,n,k,l);
              }
              else
              {
                //sum += Cijkl(i, j, m, n) * Qmnkl(m, n, k, l);
                result(i,j,k,l) += Cijkl(i,j,m,n) * Qmnkl(m,n,k,l);
              }
            }
      //  Moose::out << "\n  q"; std::cout << i + 1 << j + 1 << k + 1 << l + 1; Moose::out << " = " << _vals[i][j][k][l] << ";";
        }
  return result;
  //Moose::out << "\n Complete.";
}



Real
electrostrictiveProduct(const RankFourTensor & qijkl, unsigned int i, const RealVectorValue & v, unsigned int k, const RealVectorValue & w)
{
  // RankFourTensor qijkl;
  //Sum over (j,l) q_ijkl * v(j) * w(l) with k = _component
  Real sum = 0.0;
  for(unsigned int j = 0; j < LIBMESH_DIM; ++j)
    for(unsigned int l = 0; l < LIBMESH_DIM; ++l)
    {
      sum += qijkl(i, j, k, l) * v(j) * w(l);
    }
  return sum;
}

Real
electrostrictiveProduct(const RankFourTensor & qijkl, unsigned int i, const RealVectorValue & v, unsigned int k, const unsigned int l)
{
  // RankFourTensor qijkl;
  //Sum over j q_ijkl * v(j) where k and l = _component (used for DiagJacobian)
  Real sum = 0.0;
  for(unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      sum += qijkl(i,j,k,l) * v(j);
    }
  return sum;
}

}
