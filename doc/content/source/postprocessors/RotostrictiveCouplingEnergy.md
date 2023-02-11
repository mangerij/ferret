# RotostrictiveCouplingEnergy

!syntax description /Postprocessors/RotostrictiveCouplingEnergy

## Overview

Computes the rotostrictive coupling energy in the simulation box (volume $\Omega$),

\begin{equation}
  \begin{aligned}
    F_\mathrm{rotostr} &= \int\limits_\Omega d^3 \mathbf{r} \left\{r_{11} \left(A_x^2 \varepsilon_{xx} + A_y^2 \varepsilon_{yy} + A_z^2 \varepsilon_{zz}\right)\right\} \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left\{ r_{12} \left[ \left(A_x^2 + A_y^2\right) \varepsilon_{zz} + \left(A_y^2 + A_z^2\right) \varepsilon_{xx} + \left(A_x^2 + A_z^2\right) \varepsilon_{yy}\right]\right\} \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left\{ r_{44} \left(A_x A_y \varepsilon_{yx} + A_x A_z \varepsilon_{xz} + A_y A_z \varepsilon_{yz}\right)\right\},
  \end{aligned}
\end{equation}

suitable for cubic parent phase ferroelectric materials with antiferrodistortive modes (i.e. $\mathrm{BiFeO}_3$). The coefficients $r_{11}, r_{12}$, and $r_{44}$ are the rotostrictive tensor components (in Voight notation). The vector $\mathbf{A}$ is the oxygen octahedral antiphase tilt field and $\varepsilon_{ij}$ the elastic strain tensor of rank two. The units of the energy are given by units of $\mathbf{A}$ and $r_{ij}$ but an optional input parameter flag `energy`$\_$`scale` is provided allowing to quickly convert the postprocessed units.

## Example Input File Syntax

!! Describe and include an example of how to use the RotostrictiveCouplingEnergy object.

!syntax parameters /Postprocessors/RotostrictiveCouplingEnergy

!syntax inputs /Postprocessors/RotostrictiveCouplingEnergy

!syntax children /Postprocessors/RotostrictiveCouplingEnergy
