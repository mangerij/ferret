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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/


#pragma once

#include "TimeIntegrator.h"

class RungeKuttaFehlberg56;

template <>
InputParameters validParams<RungeKuttaFehlberg56>();


class RungeKuttaFehlberg56 : public TimeIntegrator
{
public:
  static InputParameters validParams();

  RungeKuttaFehlberg56(const InputParameters & parameters);

  virtual int order() override { return 4; }
  virtual void computeTimeDerivatives() override;
  virtual void computeADTimeDerivatives(DualReal & ad_u_dot,
                                        const dof_id_type & dof,
                                        DualReal & ad_u_dotdot) const override;
  virtual void solve() override;
  virtual void postResidual(NumericVector<Number> & residual) override;

protected:
  /**
   * Helper function that actually does the math for computing the time derivative
   */
  template <typename T, typename T2>
  void computeTimeDerivativeHelper(T & u_dot, const T2 & u_old) const;

  // Indicates the current stage.
  unsigned int _stage;

  // The number of stages in the method.  According to S9.4.2/4 of the
  // standard, we can specify a constant initializer like this for
  // integral types, it does not have to appear outside the class
  // definition.
  static const unsigned int _n_stages = 8;

  // Store pointers to the various stage residuals
  NumericVector<Number> * _stage_residuals[_n_stages];

  // Butcher tableau "C" parameters derived from _gamma
  static const Real _c[_n_stages];

  // Butcher tableau "A" values derived from _gamma.  We only use the
  // lower triangle of this.
  static const Real _a[_n_stages][_n_stages];
};

template <typename T, typename T2>
void
RungeKuttaFehlberg56::computeTimeDerivativeHelper(T & u_dot, const T2 & u_old) const
{
  u_dot -= u_old;
  u_dot *= 1. / _dt;
}
