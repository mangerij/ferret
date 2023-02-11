# AFMSublatticeSuperexchangeEnergyDensity

!syntax description /AuxKernels/AFMSublatticeSuperexchangeEnergyDensity

## Overview

Computes the free energy density due to the short-range superexchange interaction,

\begin{equation}
  \begin{aligned}
    f_\mathrm{super} = D_e \left(\mathbf{L}^2 - \mathbf{m}^2\right)
  \end{aligned}
\end{equation}

where $D_e$ is the coupling coefficient driving collinearity of the two sublattices $\mathbf{m}_1$ and $\mathbf{m}_2$. The quantities $\mathbf{L}$ and $\mathbf{m}$ are the antiferromagnetic Néel and total magnetization vectors respectively. Their definitions are $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 - \mathbf{m}_2$.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMSublatticeSuperexchangeEnergyDensity object.

!syntax parameters /AuxKernels/AFMSublatticeSuperexchangeEnergyDensity

!syntax inputs /AuxKernels/AFMSublatticeSuperexchangeEnergyDensity

!syntax children /AuxKernels/AFMSublatticeSuperexchangeEnergyDensity
