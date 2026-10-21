# MasterMagneticDMIEnergy

!syntax description /Postprocessors/MasterMagneticDMIEnergy

## Overview

Computes the free energy of the simulation box (volume $\Omega$) corresponding to the Dzyaloshinskii-Moriya
interaction of a ferromagnet,

\begin{equation}
  \begin{aligned}
    F_\mathrm{DMI}^\mathrm{bulk} &= \int\limits_\Omega d^3\mathbf{r} \, D\, \mathbf{m} \cdot \left(\nabla \times \mathbf{m}\right), \\
    F_\mathrm{DMI}^\mathrm{interfacial} &= \int\limits_\Omega d^3\mathbf{r} \, D \left[ m_z \left(\nabla \cdot \mathbf{m}\right) - \left(\mathbf{m}\cdot\nabla\right) m_z \right],
  \end{aligned}
\end{equation}

selected by `dmi_type`, with the same sign convention as the kernel [MasterDMICartLLG](MasterDMICartLLG.md). The units
of the energy are given by units of $D$ times the mesh volume, and the optional input parameter `energy`$\_$`scale`
converts them (e.g. $10^{-27}$ for a mesh in nm with `D` supplied as $D/10^{-9}$, giving J).

## Example Input File Syntax

The DMI class is taken from `dmi_type` (set in `[GlobalParams]` in this example):

!listing test/tests/magnets/sp6_1d.i block=Postprocessors/E_dmi

It enters the total energy with the other terms:

!listing test/tests/magnets/sp6_1d.i block=Postprocessors/E_total

!syntax parameters /Postprocessors/MasterMagneticDMIEnergy

!syntax inputs /Postprocessors/MasterMagneticDMIEnergy

!syntax children /Postprocessors/MasterMagneticDMIEnergy
