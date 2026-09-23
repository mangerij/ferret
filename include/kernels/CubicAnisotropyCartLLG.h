/**
 * CubicAnisotropyCartLLG -- LLG torque from CUBIC magnetocrystalline anisotropy of a single
 * (ferromagnetic) lattice.  Generated from the unified_LT free energy of the mGM4+ order
 * parameter of Pm-3m 1' at fixed |m| = 1, where the anisotropic quartic sector has exactly one
 * independent invariant:
 *
 *     F = Kc1 (mx^2 my^2 + my^2 mz^2 + mz^2 mx^2)          Kc1 > 0: <100> easy (e.g. bcc Fe)
 *                                                          Kc1 < 0: <111> easy (e.g. Ni)
 *
 * Residual (same convention as MasterAnisotropyCartLLG, verified symbolically):
 *     R_i = -(g0/(mu0 Hscale Ms)) [m x g + alpha m x (m x g)]_i psi / (1 + alpha^2),  g = dF/dm
 *
 * Material properties: alpha, Kc1, Ms.   Params: g0, mu0, Hscale.
 * The FM analogue of AFMSingleIonCubicSixthAnisotropy (which is the SIXTH-order cubic
 * invariant on two sublattices).
 */

#pragma once
#include "Kernel.h"

class CubicAnisotropyCartLLG : public Kernel
{
public:
  static InputParameters validParams();
  CubicAnisotropyCartLLG(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const unsigned int _component;
  const unsigned int _mag_x_var;
  const unsigned int _mag_y_var;
  const unsigned int _mag_z_var;
  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _Kc1;
  const MaterialProperty<Real> & _Ms;
  const Real _g0;
  const Real _mu0;
  const Real _Hscale;
};
