# CubicParentFEPhaseFieldAction

!syntax description /Ferret/CubicParentFEPhaseField/CubicParentFEPhaseFieldAction

## Overview

Sets up a cubic-parent-phase ferroelectric problem: the three polarization components and
the electrostatic potential, the Landau bulk and gradient kernels, the optional
electrostatic and electrostrictive couplings, the constant material properties and the
free energy postprocessors.

The mechanics are deliberately left to the SolidMechanics `QuasiStatic` physics, which
owns the displacement variables, the strain calculator and the stress divergence kernels.
The two actions are coupled through the eigenstrain: this action adds a
`ComputeCubicParentElectrostrictiveStrain` named by `eigenstrain_name`, and that same name
must appear in the SolidMechanics action's `eigenstrain_names`.

A prescribed mean strain (an epitaxial misfit, say) is supplied as a
`GenericConstantRankTwoTensor` and reaches the strain calculator through
`[GlobalParams] global_strain`. It contributes to `total_strain`, which is the quantity
the cubic-parent polarization kernels read.

`electrostatics` and `elastic` are independent, so the action covers the purely
polar-elastic problem as well as the full polar-elastic-electric one.

## Elastic energy

Two routes are provided, and `Ftotal` uses the second when it is available:

- `Felastic`, from `CubicParentElasticEnergy`, which evaluates the cubic-parent form
  directly and subtracts only its *own* eigenstrain.
- `Felastic_true`, an `ElementIntegralVariablePostprocessor` over the `f_el` aux variable
  that `ElasticEnergyAux` fills with 1/2 sigma:elastic_strain. Because the strain
  calculator has already removed *every* eigenstrain listed in the SolidMechanics action's
  `eigenstrain_names`, this is the true elastic energy no matter how many contribute one,
  with no double counting and no cross term left over.

With a single eigenstrain the two agree to roundoff. They diverge as soon as anything else
contributes an eigenstrain, and only `Felastic_true` stays correct; set
`add_elastic_energy_aux = false` to drop it and fall back to `Felastic`.

!alert note
At `t = 0` the `f_el` aux reads a material state in which the eigenstrain has not yet been
applied, so `Felastic_true` reports 1/2 eps:C:eps of the prescribed mean strain alone.
Every subsequent step is correct. This is MOOSE's legacy material-output behaviour on
`INITIAL`, not a property of this action.

## Example Input File Syntax

!listing test/tests/action/PTO_2D_E_action.i block=Ferret

Alongside the SolidMechanics physics:

!listing test/tests/action/PTO_2D_E_action.i block=Physics

!syntax parameters /Ferret/CubicParentFEPhaseField/CubicParentFEPhaseFieldAction
