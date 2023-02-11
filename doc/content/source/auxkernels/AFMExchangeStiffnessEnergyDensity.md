# AFMExchangeStiffnessEnergyDensity

!syntax description /AuxKernels/AFMExchangeStiffnessEnergyDensity

## Overview

Computes the free energy density due to the long-range exchange interactions in the antiferromagnet,

\begin{equation}
  \begin{aligned}
    f_\mathrm{\nabla L} = A_e \left[\left(\nabla L_x\right)^2 + \left(\nabla L_y\right)^2 \left(\nabla L_z\right)^2\right].
  \end{aligned}
\end{equation}

The quantity $\mathbf{L}$ is the antiferromagnetic Néel vector along with coupling constant $A_e$.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMExchangeStiffnessEnergyDensity object.

!syntax parameters /AuxKernels/AFMExchangeStiffnessEnergyDensity

!syntax inputs /AuxKernels/AFMExchangeStiffnessEnergyDensity

!syntax children /AuxKernels/AFMExchangeStiffnessEnergyDensity
