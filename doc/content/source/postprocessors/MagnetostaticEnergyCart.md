# MagnetostaticEnergyCart

!syntax description /Postprocessors/MagnetostaticEnergyCart

## Overview

Computes the free energy of the simulation box (volume $\Omega$) corresponding to magnetostatic energy,

\begin{equation}
  \begin{aligned}
    F_\mathrm{mag} = \frac{\mu_0 M_s}{2} \int\limits_\Omega d^3\mathbf{r} \, \left\{\mathbf{m} \cdot \nabla \Phi_\mathrm{H} \right\}
  \end{aligned}
\end{equation}

where $\mu_0$ is the permeability of vacuum, $M_s$ is the saturation magnetization density, $\mathbf{m}$ the normalized magnetization and $\Phi_\mathrm{H}$ the magnetostatic potential. An optional input parameter flag `energy`$\_$`scale` is provided allowing to quickly convert the postprocessed units.

## Example Input File Syntax

!! Describe and include an example of how to use the MagnetostaticEnergyCart object.

!syntax parameters /Postprocessors/MagnetostaticEnergyCart

!syntax inputs /Postprocessors/MagnetostaticEnergyCart

!syntax children /Postprocessors/MagnetostaticEnergyCart
