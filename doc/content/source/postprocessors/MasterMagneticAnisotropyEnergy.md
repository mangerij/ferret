# MasterMagneticAnisotropyEnergy


!syntax description /Postprocessors/MasterMagneticAnisotropyEnergy

## Overview

Computes the free energy of the simulation box (volume $\Omega$) corresponding to magnetocrystalline anisotropy,

\begin{equation}
  \begin{aligned}
    F_\mathrm{anis} = - \int\limits_\Omega d^3\mathbf{r} \,\, \left[ K_1  \, \left( \mathbf{m} \cdot \mathbf{w} \right)^2 + K_2  \, \left( \mathbf{m} \cdot \mathbf{w} \right)^4 \right]
  \end{aligned}
\end{equation}

where $\mathbf{m}$ is the normalized magnetization vector and $K_1$, $K_2$ the anisotropy constants ($K_1 > 0$: easy axis along $\mathbf{w}$, the convention of the kernel [MasterAnisotropyCartLLG](MasterAnisotropyCartLLG.md)). The material property `K2` is optional and is taken as zero when no material provides it. The vector $\mathbf{w}$ is the normalized director of the anisotropy. Typically, this is denoted as $\mathbf{n}$ in the literature but in order to avoid confusion with the commonly referred surface normal in the MOOSE code, we choose $\mathbf{w}$. The units of the energy are given by units of $K_{1}$ but an optional input parameter flag `energy`$\_$`scale` is provided allowing to quickly convert the postprocessed units.

## Example Input File Syntax

!! Describe and include an example of how to use the MasterMagneticAnisotropyEnergy object.

!syntax parameters /Postprocessors/MasterMagneticAnisotropyEnergy

!syntax inputs /Postprocessors/MasterMagneticAnisotropyEnergy

!syntax children /Postprocessors/MasterMagneticAnisotropyEnergy
