# AFMSublatticeDMInteractionEnergyDensity

!syntax description /AuxKernels/AFMSublatticeDMInteractionEnergyDensity

## Overview

Computes the free energy density due to the Dzyaloshinskii–Moriya interaction (DMI),

\begin{equation}
  \begin{aligned}
    f_\mathrm{DMI} = 4 D_0 \mathbf{A} \cdot \left( \mathbf{L} \times \mathbf{m}\right)
  \end{aligned}
\end{equation}

where $\mathbf{A}$, $\mathbf{L}$, and $\mathbf{m}$ are the oxygen octahedral antiphase tilts (see $\mathrm{BiFeO}_3$), antiferromagnetic Néel, and total magnetization vectors. The quantity $D_0$ is the DMI coupling strength in units of energy density divided by units of $\mathbf{A}$. The factor of $4$ is our convention and the quantities $\mathbf{L}$ and $\mathbf{m}$ as assumed to be normalized such that $|\mathbf{L}|+|\mathbf{m}| = 2$ with $|\mathbf{L}| \gg |\mathbf{m}|$. By setting $\mathbf{A}$ to some other vector relevant to the problem, this can be generalized to other antiferromagnets.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMSublatticeDMInteractionEnergyDensity object.

!syntax parameters /AuxKernels/AFMSublatticeDMInteractionEnergyDensity

!syntax inputs /AuxKernels/AFMSublatticeDMInteractionEnergyDensity

!syntax children /AuxKernels/AFMSublatticeDMInteractionEnergyDensity
