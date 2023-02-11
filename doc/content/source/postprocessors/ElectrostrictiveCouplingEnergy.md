# ElectrostrictiveCouplingEnergy

!syntax description /Postprocessors/ElectrostrictiveCouplingEnergy

## Overview

Computes the electrostrictive coupling energy in the simulation box (volume $\Omega$),

\begin{equation}
  \begin{aligned}
    F_\mathrm{elecstr} &= \int\limits_\Omega d^3 \mathbf{r} \left\{q_{11} \left(P_x^2 \varepsilon_{xx} + P_y^2 \varepsilon_{yy} + P_z^2 \varepsilon_{zz}\right)\right\} \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left\{ q_{12} \left[ \left(P_x^2 + P_y^2\right) \varepsilon_{zz} + \left(P_y^2 + P_z^2\right) \varepsilon_{xx} + \left(P_x^2 + P_z^2\right) \varepsilon_{yy}\right]\right\} \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left\{ q_{44} \left(P_x P_y \varepsilon_{yx} + P_x P_z \varepsilon_{xz} + P_y P_z \varepsilon_{yz}\right)\right\},
  \end{aligned}
\end{equation}

suitable for cubic parent phase ferroelectric materials (i.e. $\mathrm{PbTiO}_3$ and $\mathrm{BaTiO}_3$). The coefficients $q_{11}, q_{12}$, and $q_{44}$ are the electrostrictive tensor components (in Voight notation). The vector $\mathbf{P}$ is the polarization field and $\varepsilon_{ij}$ the elastic strain tensor of rank two. The units of the energy are given by units of $\mathbf{P}$ and $q_{ij}$ but an optional input parameter flag `energy`$\_$`scale` is provided allowing to quickly convert the postprocessed units.


## Example Input File Syntax

!! Describe and include an example of how to use the ElectrostrictiveCouplingEnergy object.

!syntax parameters /Postprocessors/ElectrostrictiveCouplingEnergy

!syntax inputs /Postprocessors/ElectrostrictiveCouplingEnergy

!syntax children /Postprocessors/ElectrostrictiveCouplingEnergy
