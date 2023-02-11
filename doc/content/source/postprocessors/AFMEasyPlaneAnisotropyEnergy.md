# AFMEasyPlaneAnisotropyEnergy

!syntax description /Postprocessors/AFMEasyPlaneAnisotropyEnergy

## Overview

Computes the free energy of the simulation box (volume $\Omega$) corresponding to easy-plane magnetic anisotropy (in the case of $\mathrm{BiFeO}_3$),

\begin{equation}
  \begin{aligned}
    F_\mathrm{easy} = \int\limits_\Omega d^3\mathbf{r}  \left[\sum\limits_{\eta = 1}^2 K_1 \left(\textbf{m}_\eta \cdot \hat{\mathbf{P}}\right)^2\right]
  \end{aligned}
\end{equation}

where $\mathbf{m}_\eta$ are the magnetic sublattices $\eta = 1,2$ and $K_1 < 0$ the anisotropy constant. The vector $\hat{\mathbf{P}}$ is the normalized director of the ferroelectric polarization.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMEasyPlaneAnisotropyEnergy object.

!syntax parameters /Postprocessors/AFMEasyPlaneAnisotropyEnergy

!syntax inputs /Postprocessors/AFMEasyPlaneAnisotropyEnergy

!syntax children /Postprocessors/AFMEasyPlaneAnisotropyEnergy
