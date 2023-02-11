# AFMEasyPlaneAnisotropyEnergyDensity

!syntax description /AuxKernels/AFMEasyPlaneAnisotropyEnergyDensity

## Overview

Computes the free energy density due to magnetic anisotropy in the ferroelectric antiferromagnet $\mathrm{BiFeO}_3$ corresponding to

\begin{equation}
  \begin{aligned}
    f_\mathrm{easy} = \sum\limits_{\eta = 1}^2 K_1 \left(\textbf{m}_\eta \cdot \hat{\mathbf{P}}\right)^2
  \end{aligned}
\end{equation}

where $\mathbf{m}_\eta$ are the magnetic sublattices $\eta = 1,2$ and $K_1 < 0$ the anisotropy constant. The vector $\hat{\mathbf{P}}$ is the normalized director of the ferroelectric polarization.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMEasyPlaneAnisotropyEnergyDensity object.

!syntax parameters /AuxKernels/AFMEasyPlaneAnisotropyEnergyDensity

!syntax inputs /AuxKernels/AFMEasyPlaneAnisotropyEnergyDensity

!syntax children /AuxKernels/AFMEasyPlaneAnisotropyEnergyDensity
