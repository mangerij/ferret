# MasterMagneticZeemanEnergyCart

!syntax description /Postprocessors/MasterMagneticZeemanEnergyCart

## Overview

Computes the free energy in the simulation box (volume $\Omega$) corresponding to the Zeeman interaction,

\begin{equation}
  \begin{aligned}
    F_\mathrm{Zeeman} &= - \mu_0 M_s \int\limits_\Omega d^3 \mathbf{r} \left\{ \mathbf{m}\cdot\mathbf{H}_\mathrm{ext}\right\}. \\
  \end{aligned}
\end{equation}

where $\mathbf{m}$ is the normalized magnetization vector, $M_s$ is the saturation magnetization density, $\mu_0$ is the permeability of vacuum. The applied (external) field is $\mathbf{H}_\mathrm{ext}$. The units of the energy are given by units of $\mathbf{H}_\mathrm{ext}$, $\mu_0$, and $M_s$ but an optional input parameter flag `energy`$\_$`scale` is provided allowing to quickly convert the postprocessed units.

## Example Input File Syntax

!! Describe and include an example of how to use the MasterMagneticZeemanEnergyCart object.

!syntax parameters /Postprocessors/MasterMagneticZeemanEnergyCart

!syntax inputs /Postprocessors/MasterMagneticZeemanEnergyCart

!syntax children /Postprocessors/MasterMagneticZeemanEnergyCart
