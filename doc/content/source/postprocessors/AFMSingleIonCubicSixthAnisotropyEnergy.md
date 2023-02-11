# AFMSingleIonCubicSixthAnisotropyEnergy

!syntax description /Postprocessors/AFMSingleIonCubicSixthAnisotropyEnergy

## Overview

Computes the free energy in the simulation box (volume $\Omega$) due to weak sixth order magnetic anisotropy in the ferroelectric antiferromagnet $\mathrm{BiFeO}_3$ corresponding to

\begin{equation}
  \begin{aligned}
    F_\mathrm{anis} &= \int\limits_\Omega d^3 \mathbf{r} \left[ \sum\limits_{\eta = 1}^2 \left(K_1^c + a|\mathbf{A}|^2\right)\left(m_{\eta,x}^2 m_{\eta,y}^2 m_{\eta,z}^2\right)\right]
  \end{aligned}
\end{equation}

where $\mathbf{m}_\eta$ are the magnetic sublattices $\eta = 1,2$ and $K_1^c$ and $a$ coupling constants.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMSingleIonCubicSixthAnisotropyEnergy object.

!syntax parameters /Postprocessors/AFMSingleIonCubicSixthAnisotropyEnergy

!syntax inputs /Postprocessors/AFMSingleIonCubicSixthAnisotropyEnergy

!syntax children /Postprocessors/AFMSingleIonCubicSixthAnisotropyEnergy
