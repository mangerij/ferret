//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RankSixTensor.h"

// MOOSE includes
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MooseEnum.h"
#include "MooseException.h"
#include "MooseUtils.h"
#include "MatrixTools.h"
#include "MaterialProperty.h"
#include "PermutationTensor.h"

#include "libmesh/utility.h"

// C++ includes
#include <iomanip>
#include <ostream>

template <>
void
mooseSetToZero<RankSixTensor>(RankSixTensor & v)
{
  v.zero();
}

template <>
void
dataStore(std::ostream & stream, RankSixTensor & rft, void * context)
{
  dataStore(stream, rft._vals, context);
}

template <>
void
dataLoad(std::istream & stream, RankSixTensor & rft, void * context)
{
  dataLoad(stream, rft._vals, context);
}

MooseEnum
RankSixTensor::fillMethodEnum()
{
  return MooseEnum("general");
}

RankSixTensor::RankSixTensor()
{
  mooseAssert(N == 3, "RankSixTensor is currently only tested for 3 dimensions.");

  unsigned int index = 0;
  for (unsigned int i = 0; i < N6; ++i)
    _vals[index++] = 0.0;
}

RankSixTensor::RankSixTensor(const InitMethod init)
{
  unsigned int index = 0;
  switch (init)
  {
    case initNone:
      break;

    case initIdentity:
      zero();
      for (unsigned int i = 0; i < N; ++i)
        (*this)(i, i, i, i, i, i) = 1.0;
      break;

    case initIdentitySix:
      for (unsigned int i = 0; i < N; ++i)
        for (unsigned int j = 0; j < N; ++j)
          for (unsigned int k = 0; k < N; ++k)
            for (unsigned int l = 0; l < N; ++l)
              for (unsigned int m = 0; m < N; ++m)
                for (unsigned int n = 0; n < N; ++n)
              _vals[index++] = (i == k) && (j == l) && (m == n);
      break;

    case initIdentitySymmetricSix:
      for (unsigned int i = 0; i < N; ++i)
        for (unsigned int j = 0; j < N; ++j)
          for (unsigned int k = 0; k < N; ++k)
            for (unsigned int l = 0; l < N; ++l)
              for (unsigned int m = 0; m < N; ++m)
                for (unsigned int n = 0; n < N; ++n)
              _vals[index++] = 0.5 * ((i == k) && (j == l) && (m == n)) + 0.5 * ((i == l) && (j == k) && (m == n));
      break;

    default:
      mooseError("Unknown RankSixTensor initialization pattern.");
  }
}

RankSixTensor::RankSixTensor(const std::vector<Real> & input, FillMethod fill_method)
{
  fillFromInputVector(input, fill_method);
}

void
RankSixTensor::zero()
{
  for (unsigned int i = 0; i < N6; ++i)
    _vals[i] = 0.0;
}

RankSixTensor &
RankSixTensor::operator=(const RankSixTensor & a)
{
  for (unsigned int i = 0; i < N6; ++i)
    _vals[i] = a._vals[i];
  return *this;
}

//RankFourTensor RankSixTensor::operator*(const RankFourTensor & b) const
//{
//  RankFourTensor result;
//
//  unsigned int index = 0;
//  for (unsigned int ij = 0; ij < N4; ++ij)
//  {
//    Real tmp = 0;
//    for (unsigned int kl = 0; kl < N4; ++kl)
//      tmp += _vals[index++] * b._coords[kl];
//    result._coords[ij] = tmp;
//  }
//
//  return result;
//}

RankSixTensor RankSixTensor::operator*(const Real b) const
{
  RankSixTensor result;

  for (unsigned int i = 0; i < N6; ++i)
    result._vals[i] = _vals[i] * b;

  return result;
}

RankSixTensor &
RankSixTensor::operator*=(const Real a)
{
  for (unsigned int i = 0; i < N6; ++i)
    _vals[i] *= a;
  return *this;
}

RankSixTensor
RankSixTensor::operator/(const Real b) const
{
  RankSixTensor result;
  for (unsigned int i = 0; i < N6; ++i)
    result._vals[i] = _vals[i] / b;
  return result;
}

RankSixTensor &
RankSixTensor::operator/=(const Real a)
{
  for (unsigned int i = 0; i < N6; ++i)
    _vals[i] /= a;
  return *this;
}

RankSixTensor &
RankSixTensor::operator+=(const RankSixTensor & a)
{
  for (unsigned int i = 0; i < N6; ++i)
    _vals[i] += a._vals[i];
  return *this;
}

RankSixTensor
RankSixTensor::operator+(const RankSixTensor & b) const
{
  RankSixTensor result;
  for (unsigned int i = 0; i < N6; ++i)
    result._vals[i] = _vals[i] + b._vals[i];
  return result;
}

RankSixTensor &
RankSixTensor::operator-=(const RankSixTensor & a)
{
  for (unsigned int i = 0; i < N6; ++i)
    _vals[i] -= a._vals[i];
  return *this;
}

RankSixTensor
RankSixTensor::operator-(const RankSixTensor & b) const
{
  RankSixTensor result;
  for (unsigned int i = 0; i < N6; ++i)
    result._vals[i] = _vals[i] - b._vals[i];
  return result;
}

RankSixTensor
RankSixTensor::operator-() const
{
  RankSixTensor result;
  for (unsigned int i = 0; i < N6; ++i)
    result._vals[i] = -_vals[i];
  return result;
}

Real
RankSixTensor::L2norm() const
{
  Real l2 = 0;

  for (unsigned int i = 0; i < N6; ++i)
    l2 += Utility::pow<2>(_vals[i]);

  return std::sqrt(l2);
}

void
RankSixTensor::rotate(const RealTensorValue & R)
{
  RankSixTensor old = *this;

  unsigned int index = 0;
  for (unsigned int i = 0; i < N; ++i)
  {
    for (unsigned int j = 0; j < N; ++j)
    {
      for (unsigned int k = 0; k < N; ++k)
      {
        for (unsigned int l = 0; l < N; ++l)
        {
          for (unsigned int m = 0; m < N; ++m)
          {
            for (unsigned int n = 0; m < N; ++n)
            {
              unsigned int index2 = 0;
              Real sum = 0.0;
              for (unsigned int r = 0; r < N; ++r)
                {
                  const Real a = R(i, r);
                  for (unsigned int s = 0; s < N; ++s)
                  {
                    const Real ab = a * R(j, s);
                    for (unsigned int o = 0; o < N; ++o)
                    {
                      const Real abc = ab * R(k, o);
                      for (unsigned int p = 0; p < N; ++p)
                      {
                        const Real abcd = abc * R(j, p);
                        for (unsigned int q = 0; q < N; ++q)
                        {
                          const Real abcde = abcd * R(k, q);
                          for (unsigned int z = 0; z < N; ++z)
                            sum += abcde * R(l, z) * old._vals[index2++];
                        }
                      }
                    }
                  }
                }
            _vals[index++] = sum;
          }
        }
      }
    }
  }
}
}

void
RankSixTensor::fillFromInputVector(const std::vector<Real> & input, FillMethod fill_method)
{
  zero();

  switch (fill_method)
  {
    case general:
      fillGeneralFromInputVector(input);
      break;
    default:
      mooseError("fillFromInputVector called with unknown fill_method of ", fill_method);
  }
}

void
RankSixTensor::fillGeneralFromInputVector(const std::vector<Real> & input)
{
  if (input.size() != 729)
    mooseError("To use fillGeneralFromInputVector, your input must have size 729. Yours has size ",
               input.size());

  for (unsigned int i = 0; i < N6; ++i)
    _vals[i] = input[i];
}
