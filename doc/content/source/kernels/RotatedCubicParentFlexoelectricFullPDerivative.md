# RotatedCubicParentFlexoelectricFullPDerivative

!syntax description /Kernels/RotatedCubicParentFlexoelectricFullPDerivative

## Overview

The elastic energy density of a cubic parent phase with electrostrictive and flexoelectric
eigenstrains is

\begin{equation}
  f_\mathrm{elastic} = \frac{1}{2} C_{ijkl}\left(\varepsilon_{ij} - \varepsilon^Q_{ij} - \varepsilon^F_{ij}\right)\left(\varepsilon_{kl} - \varepsilon^Q_{kl} - \varepsilon^F_{kl}\right),
\end{equation}

with the electrostrictive and flexoelectric eigenstrains

\begin{equation}
  \varepsilon^Q_{ij} = Q_{ijkl} P_k P_l, \qquad \varepsilon^F_{ij} = -A_{ijkl} \frac{\partial P_k}{\partial x_l}, \qquad A_{ijkl} = s_{ijmn} f_{mnkl},
\end{equation}

where $\varepsilon_{ij}$ is the total strain, $s_{ijkl}$ the elastic compliance and $f_{ijkl}$ the
flexocoupling tensor. This form is the exact Legendre transform of the Gibbs functional of
[!cite](Gu2014), including the $\frac{1}{2} C \varepsilon^F \varepsilon^F$ term.

This kernel computes the residual of the polarization equation due to this energy,

\begin{equation}
  \mathcal{R}_{P_i} = \left(\psi, \frac{\partial f_\mathrm{elastic}}{\partial P_i}\right) + \left(\frac{\partial \psi}{\partial x_l}, \frac{\partial f_\mathrm{elastic}}{\partial (\partial P_i/\partial x_l)}\right),
\end{equation}

which contains +both+ the electrostrictive and the flexoelectric contributions. It is used
+instead of+ [CubicParentElasticPDerivative](CubicParentElasticPDerivative.md), never together
with it; with $f = 0$ the two are identical.

The residual, the diagonal and off-diagonal $P$ Jacobian, the $P \leftarrow u$ block (when
`displacements` is given) and the $P \leftrightarrow \bar{\varepsilon}$ blocks with a SCALAR
mean strain (when `global_strain` is given) are all generated from the single energy density
above, so they are exact derivatives of each other and of the energy computed by the
corresponding postprocessor.

!alert note title=Conventions
`mu11`, `mu12` and `mu44` are the +tensor+ components $f_{1111}$, $f_{1122}$ and $f_{1212}$ of the
flexocoupling tensor; the compliances `s11`, `s12`, `s44` are in Voigt form ($s_{1212} = s_{44}/4$).
With the sign above, the flexocoupling $f = C:F$ computed from the published $F$ of
[!cite](Gu2014) is entered directly and gives their tail-to-tail Néel component. The flexoelectric
eigenstrain is +added+ to the electrostrictive one (no `CompositeEigenstrain` weight of $-1$).
The deprecated NV/RI flexoelectric objects of earlier development versions carried the shear
slot as $2 s_{44}\mu_{44}$, i.e. four times too large, and expected the opposite sign.

!alert warning title=Literature flexocoupling values
Published flexocoupling tensors come in several index conventions. In particular the
displacement-gradient form $P_i = \mu_{ijkl} u_{j,kl}$ of [!cite](Maranganti2009) mixes the 12
and 44 components relative to the strain-gradient form used here, and Voigt shear factors of
1, 2 or 4 appear in the literature. Convert to the tensor components $f_{1111}, f_{1122}, f_{1212}$
of the strain-gradient form before entering them.

!alert note title=Rotated crystal axes
The crystal axes are turned about $z$ by the material property `twist` (radians); the polarization,
its gradient and the strain are expressed in the crystal frame, the cubic expressions are evaluated
there and the result is turned back. For the displacement equation to see the same crystal,
`ComputeElasticityTensor` must be given `euler_angle_3` $= -$`twist` (in degrees).


## Example Input File Syntax

!listing test/tests/flexo/flexo_180DW_theta15.i block=Kernels/flexo_polar_coupled_x

!syntax parameters /Kernels/RotatedCubicParentFlexoelectricFullPDerivative

!syntax inputs /Kernels/RotatedCubicParentFlexoelectricFullPDerivative

!syntax children /Kernels/RotatedCubicParentFlexoelectricFullPDerivative
