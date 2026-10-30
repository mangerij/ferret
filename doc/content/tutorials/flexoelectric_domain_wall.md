# Flexoelectricity at a 180 degree domain wall

A 180 degree ferroelectric domain wall is, at the level of the bulk Landau energy, a pure Ising
wall: the polarization reverses along one axis and the transverse components vanish everywhere.
Gradient terms and flexoelectric coupling break that picture. Because the wall carries a large
polarization gradient, the flexoelectric term converts that gradient into a strain, and the
elastic response feeds back a *transverse* polarization. The wall acquires Bloch-like and
Neel-like character.

This tutorial follows
[`test/tests/flexo/flexo_180DW_theta15.i`](https://github.com/mangerij/ferret/blob/devel/test/tests/flexo/flexo_180DW_theta15.i),
which relaxes a 180 degree wall in BaTiO$_3$ with flexoelectric coupling switched on.

## Geometry and seed

The wall normal is $x$, and the problem is quasi one-dimensional: a long bar with a single
element across $y$ and $z$.

!listing test/tests/flexo/flexo_180DW_theta15.i
         block=Mesh
         link=False
         language=python

The initial condition seeds the wall with a cosine profile in $P_z$,

!listing test/tests/flexo/flexo_180DW_theta15.i
         block=Functions
         link=False
         language=python

with `period` chosen as $2\pi / L_x$ so that one full period spans the bar. This produces two
180 degree walls in the periodic cell, which is what the profile plot below shows.

The parameter `thetadeg` at the top of the input rotates the crystal axes relative to the wall
normal. This matters a great deal, as the results show.

## The flexoelectric coupling

Three objects carry the physics. The material computes the flexoelectric contribution to the
strain,

!listing test/tests/flexo/flexo_180DW_theta15.i
         block=Materials/flexoelectric_eigenstrain
         link=False
         language=python

the kernel supplies the corresponding $\partial / \partial P$ driving force in the polarization
equations, and a postprocessor reports the energy. All three are the `Rotated` variants, which
apply the Euler rotation set by `thetadeg` to the flexoelectric tensor.

The transverse polarization components are monitored directly:

!listing test/tests/flexo/flexo_180DW_theta15.i
         block=Postprocessors
         link=False
         language=python

`maxPy` and `minPy` are the quantity of interest - the Bloch component of the wall.

## What the relaxed wall looks like

!media media/flexo_180DW_profiles.png style=display:block;margin-left:auto;margin-right:auto;width:95%; caption=Polarization profiles across the relaxed 180 degree wall. Bottom: $P_z$ reverses across each wall. Middle: $P_y$, the Bloch component, is antisymmetric about the wall centre and is non-zero only for the rotated wall. Top: $P_x$, the Neel component. Solid curves are $\theta = 0$, dashed are $\theta = 20.2$ degrees. NV and PV denote two forms of the Voigt contraction in the flexoelectric coupling. id=fig-flexo-profiles

The structure of [fig-flexo-profiles] is the point of the tutorial:

- $P_z$ (bottom) is the ordinary 180 degree reversal, saturating at $\pm P_s \approx 0.26$
  C/m$^2$ and essentially independent of $\theta$.
- $P_y$ (middle) is the +Bloch component+. It is antisymmetric about the wall centre, peaks at
  roughly $3 \times 10^{-3}$ C/m$^2$, and is *identically zero for the unrotated wall*. It exists
  only because the rotation makes the flexoelectric tensor couple the $P_z$ gradient to a
  transverse strain.
- $P_x$ (top) is the Neel component and behaves similarly, though it is about a factor of three
  smaller here.

So the flexoelectric coupling does not merely perturb the wall - for $\theta \neq 0$ it changes
its character, from Ising to mixed Ising-Bloch-Neel.

!alert note title=The figure comes from a longer run than the regression test
`flexo_180DW_theta15.i` is a regression test and stops after a few steps, long before the
transverse components have grown to their relaxed values - its `maxPy` only reaches about
$5 \times 10^{-5}$. The profiles above come from fully converged runs of the same model, in a
larger cell, at $\theta = 0$ and $\theta = 20.2$ degrees, where the `Terminator` on
`perc_change` has fired. To reproduce them, raise `num_steps` and `end_time` and let the energy
rate converge.

## A caution on conventions

The two colours in [fig-flexo-profiles] are two different ways of applying the Voigt contraction
to the flexoelectric tensor. They give the same $P_z$ profile but transverse components that
differ by tens of percent. The shear factor in the Voigt convention is easy to get wrong and the
error is invisible in the Ising component - it shows up only in exactly the small transverse
signal this tutorial is about. If you implement a new flexoelectric coupling, check it against a
rotated wall rather than an unrotated one.

!content pagination previous=tutorials/spin_waves.md next=tutorials/thermoelectric.md
