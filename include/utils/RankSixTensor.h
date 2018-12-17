//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef RANKSIXTENSOR_H
#define RANKSIXTENSOR_H

// MOOSE includes
#include "DataIO.h"

#include "libmesh/tensor_value.h"
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h"

// Forward declarations
class MooseEnum;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
template <typename>
class RankFourTensorTempl;
typedef RankFourTensorTempl<Real> RankFourTensor;
class RankSixTensor;

template <typename T>
void mooseSetToZero(T & v);

/**
 * Helper function template specialization to set an object to zero.
 * Needed by DerivativeMaterialInterface
 */
template <>
void mooseSetToZero<RankSixTensor>(RankSixTensor & v);

/**
 * RankSixTensor is designed to handle any N-dimensional sixth order tensor, C.
 *
 *
 * Since N is hard-coded to 3, RankSixTensor holds 729 separate C_ijklmn entries.
 * Within the code i = 0, 1, 2, but this object provides methods to extract the entries
 * with i = 1, 2, 3, and some of the documentation is also written in this way.
 */
class RankSixTensor
{
public:
  /// Initialization method
  enum InitMethod
  {
    initNone,
    initIdentity,
    initIdentitySix,
    initIdentitySymmetricSix
  };

  /**
   * To fill up the 81 entries in the 4th-order tensor, fillFromInputVector
   * is called with one of the following fill_methods.
   * See the fill*FromInputVector functions for more details
   */
  enum FillMethod
  {
    antisymmetric,
    symmetric9,
    symmetric21,
    general_isotropic,
    symmetric_isotropic,
    symmetric_isotropic_E_nu,
    antisymmetric_isotropic,
    axisymmetric_rz,
    general,
    principal
  };

  /// Default constructor; fills to zero
  RankSixTensor();

  /// Select specific initialization pattern
  RankSixTensor(const InitMethod);

  /// Fill from vector
  RankSixTensor(const std::vector<Real> &, FillMethod);

  // Named constructors
  static RankSixTensor Identity() { return RankSixTensor(initIdentity); }
  static RankSixTensor IdentitySix() { return RankSixTensor(initIdentitySix); };

  /// Gets the value for the index specified.  Takes index = 0,1,2
  inline Real & operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int m, unsigned int n)
  {
    return _vals[((((i * LIBMESH_DIM + j) * LIBMESH_DIM + k) * LIBMESH_DIM + l )* LIBMESH_DIM + m )* LIBMESH_DIM + n];
  }

  /**
   * Gets the value for the index specified.  Takes index = 0,1,2
   * used for const
   */
  inline Real operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int m, unsigned int n) const
  {
    return _vals[((((i * LIBMESH_DIM + j) * LIBMESH_DIM + k) * LIBMESH_DIM + l ) * LIBMESH_DIM + m )* LIBMESH_DIM + n];
  }

  /// Zeros out the tensor.
  void zero();

  /// Print the rank six tensor
  void print(std::ostream & stm = Moose::out) const;

  /// copies values from a into this tensor
  RankSixTensor & operator=(const RankSixTensor & a);

  /// C_ijklmn*a_mn
  RankFourTensor operator*(const RankFourTensor & a) const;

  /// C_ijklmn*a_klmn
  RankTwoTensor operator*(const RankTwoTensor & a) const;

  /// C_ijklmn*a
  RankSixTensor operator*(const Real a) const;

  /// C_ijklmn *= a
  RankSixTensor & operator*=(const Real a);

  /// C_ijkl/a
  RankSixTensor operator/(const Real a) const;

  /// C_ijkl /= a  for all i, j, k, l, m, n
  RankSixTensor & operator/=(const Real a);

  /// C_ijkl += a_ijkl  for all i, j, k, l, m, n
  RankSixTensor & operator+=(const RankSixTensor & a);

  /// C_ijkl + a_ijkl
  RankSixTensor operator+(const RankSixTensor & a) const;

  /// C_ijkl -= a_ijkl
  RankSixTensor & operator-=(const RankSixTensor & a);

  /// C_ijkl - a_ijkl
  RankSixTensor operator-(const RankSixTensor & a) const;

  /// -C_ijkl
  RankSixTensor operator-() const;

  /// C_ijpqmn*a_pqklrs
  //RankSixTensor operator*(const RankSixTensor & a) const;

  /// C_ijpqmn*a_pqklmn ???
  //RankSixTensor operator*(const RankSixTensor & a) const;

  /// sqrt(C_ijklmn*C_ijklmn)
  Real L2norm() const;

  /**
   * FIXME: This returns A_ijklmn such that C_ijklmn*A_klmnrs = 0.5*(de_im de_jn + de_in de_jm)
   * This routine assumes that C_ijkl = C_jikl = C_ijlk
   */
  RankSixTensor invSymm() const;

  /**
   * Rotate the tensor using
   * C_ijklrs = R_im R_in R_ko R_lp R_rq R_sf C_mnopqf
   */
  template <class T>
  void rotate(const T & R);

  /**
   * Rotate the tensor using
   * C_ijklrs = R_im R_jn R_ko R_lp R_rq R_sf C_mnopqf
   */
  void rotate(const RealTensorValue & R);

 // FIXME: /**
 //  * Transpose the tensor by swapping the first pair with the second pair of indices
 //  * @return C_klji
 //  */
 // RankSixTensor transposeMajor() const;

  /**
   * FIXME: Fills the tensor entries ignoring the last dimension (ie, C_ijkl=0 if any of i, j, k, or l =
   * 3).
   * Fill method depends on size of input
   * Input size = 2.  Then C_1111 = C_2222 = input[0], and C_1122 = input[1], and C_1212 = (input[0]
   * - input[1])/2,
   *                  and C_ijkl = C_jikl = C_ijlk = C_klij, and C_1211 = C_1222 = 0.
   * Input size = 9.  Then C_1111 = input[0], C_1112 = input[1], C_1122 = input[3],
   *                       C_1212 = input[4], C_1222 = input[5], C_1211 = input[6]
   *                       C_2211 = input[7], C_2212 = input[8], C_2222 = input[9]
   *                       and C_ijkl = C_jikl = C_ijlk
   */
 // void surfaceFillFromInputVector(const std::vector<Real> & input);

  /// Static method for use in validParams for getting the "fill_method"
  static MooseEnum fillMethodEnum();

  /**
   * fillFromInputVector takes some number of inputs to fill
   * the Rank-4 tensor.
   * @param input the numbers that will be placed in the tensor
   * @param fill_method this can be:
   *             antisymmetric (use fillAntisymmetricFromInputVector)
   *             symmetric9 (use fillSymmetricFromInputVector with all=false)
   *             symmetric21 (use fillSymmetricFromInputVector with all=true)
   *             general_isotropic (use fillGeneralIsotropicFrominputVector)
   *             symmetric_isotropic (use fillSymmetricIsotropicFromInputVector)
   *             antisymmetric_isotropic (use fillAntisymmetricIsotropicFromInputVector)
   *             axisymmetric_rz (use fillAxisymmetricRZFromInputVector)
   *             general (use fillGeneralFromInputVector)
   *             principal (use fillPrincipalFromInputVector)
   */
  void fillFromInputVector(const std::vector<Real> & input, FillMethod fill_method);

  ///@{ Vector-less fill API functions. See docs of the corresponding ...FromInputVector methods
  ///@}

protected:
  /// Dimensionality of rank-four tensor
  static constexpr unsigned int N = LIBMESH_DIM;
  static constexpr unsigned int N2 = N * N;
  static constexpr unsigned int N3 = N * N * N;
  static constexpr unsigned int N4 = N * N * N * N;
  static constexpr unsigned int N6 = N * N * N * N * N * N;

  /// The values of the rank-six tensor stored by
  /// index=((((i * LIBMESH_DIM + j) * LIBMESH_DIM + k) * LIBMESH_DIM + l) * LIBMESH_DIM + m ) * LIBMESH_DIM + n)
  Real _vals[N6];

  /**
   * fillGeneralFromInputVector takes 729 inputs to fill the Rank-4 tensor
   * No symmetries are explicitly maintained
   * @param input  C(i,j,k,l,m,n) = input[i*N*N*N*N*N + j*N*N*N*N + k*N*N*N + l*N*N + m*N + n]
   */
  void fillGeneralFromInputVector(const std::vector<Real> & input);

  template <class T>
  friend void dataStore(std::ostream &, T &, void *);

  template <class T>
  friend void dataLoad(std::istream &, T &, void *);

  friend class RankTwoTensorTempl<Real>;
  friend class RankThreeTensor;
  friend class RankFourTensorTempl<Real>;
};

template <>
void dataStore(std::ostream &, RankSixTensor &, void *);

template <>
void dataLoad(std::istream &, RankSixTensor &, void *);

inline RankSixTensor operator*(Real a, const RankSixTensor & b) { return b * a; }

template <class T>
void
RankSixTensor::rotate(const T & R)
{
  RankSixTensor old = *this;

  unsigned int index = 0;
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      for (unsigned int k = 0; k < N; ++k)
        for (unsigned int l = 0; l < N; ++l)
          for (unsigned int m = 0; m < N; ++m)
            for (unsigned int n = 0; n < N; ++n)
            {
              Real sum = 0.0;
              int index2 = 0;
              for (unsigned int r = 0; r < N; ++r)
              {
                Real a = R(i, r);
                for (unsigned int s = 0; s < N; ++s)
                {
                  Real ab = a * R(j, s);
                  for (unsigned int o = 0; o < N; ++o)
                  {
                    Real abc = ab * R(k, o);
                    for (unsigned int p = 0; p < N; ++p)
                    {
                      Real abcd = abc * R(k, p);
                      for (unsigned int q = 0; q < N; ++q)
                      {
                        Real abcde = abcd * R(k, q);
                        for (unsigned int q = 0; q < N; ++q)
                        {
                          for (unsigned int z = 0; z < N; ++z)
                          sum += abcde * R(l, z) * old._vals[index++];
                        }
                      }
                    }
                  }
                }
              }
             _vals[index++] = sum;
            }
}

#endif // RANKSIXTENSOR_H
