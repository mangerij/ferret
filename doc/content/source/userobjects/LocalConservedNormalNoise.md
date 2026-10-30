# LocalConservedNormalNoise

!syntax description /UserObjects/LocalConservedNormalNoise

## Overview

Provides a normally distributed random number (Box-Muller) at every quadrature point, regenerated at the start of each time step, for noise kernels such as `LocalConservedLangevinNoise` and the thermal-field route of [MasterInteractionCartLLGHConst](MasterInteractionCartLLGHConst.md).

- `conserved = true` (default) subtracts the volume average, so the noise integrates to zero, as needed for conserved (concentration-type) fields. Use `conserved = false` for non-conserved order parameters and for macrospins, where the uniform mode is physical.
- `white_noise = true` divides each value by $\sqrt{\mathrm{JxW}\,\Delta t}$. The noise is then space-time white with unit strength, so a source $-\psi\sqrt{2D}\,\xi$ has the fluctuation-dissipation covariance $2D\,M/\Delta t$ with the consistent mass matrix $M$. The stochastic Landau-Lifshitz-Bloch equation uses this ([Finite-temperature micromagnetics](theory/intro_LLB.md)).

MOOSE reseeds random objects at every step with `seed + timestep`. Give independent noise objects seeds that differ by more than the number of time steps; otherwise their streams repeat each other, shifted in time.

## Example Input File Syntax

!listing test/tests/magnets/sllb_ensemble.i block=UserObjects/hx_noise

!syntax parameters /UserObjects/LocalConservedNormalNoise

!syntax inputs /UserObjects/LocalConservedNormalNoise

!syntax children /UserObjects/LocalConservedNormalNoise
