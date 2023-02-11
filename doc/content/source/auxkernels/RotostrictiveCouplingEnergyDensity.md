# RotostrictiveCouplingEnergyDensity

!syntax description /AuxKernels/RotostrictiveCouplingEnergyDensity

## Overview

Computes the rotostrictive coupling energy density,

\begin{equation}
  \begin{aligned}
    f_\mathrm{rotostr} &= r_{11} \left(A_x^2 \varepsilon_{xx} + A_y^2 \varepsilon_{yy} + A_z^2 \varepsilon_{zz}\right) \\
    &+ r_{12} \left[ \left(A_x^2 + A_y^2\right) \varepsilon_{zz} + \left(A_y^2 + A_z^2\right) \varepsilon_{xx} + \left(A_x^2 + A_z^2\right) \varepsilon_{yy}\right] \\
    &+ r_{44} \left(A_x A_y \varepsilon_{yx} + A_x A_z \varepsilon_{xz} + A_y A_z \varepsilon_{yz}\right),
  \end{aligned}
\end{equation}

suitable for cubic parent phase ferroelectric materials with coupling to antiferrodistortive modes (i.e. $\mathrm{BiFeO}_3$). The coefficients $r_{11}, r_{12}$, and $r_{44}$ are the rotostrictive tensor components (in Voight notation). The vector $\mathbf{A}$ is the oxygen octahedral antiphase tilt field and $\varepsilon_{ij}$ the elastic strain tensor of rank two. The units of the energy are given by units of $\mathbf{A}$ and $r_{ij}$ but an optional input parameter flag `energy`$\_$`scale` is provided allowing to quickly convert the postprocessed units.

## Example Input File Syntax

!! Describe and include an example of how to use the RotostrictiveCouplingEnergyDensity object.

!syntax parameters /AuxKernels/RotostrictiveCouplingEnergyDensity

!syntax inputs /AuxKernels/RotostrictiveCouplingEnergyDensity

!syntax children /AuxKernels/RotostrictiveCouplingEnergyDensity
