# AFMSublatticeSuperexchangeEnergy

!syntax description /Postprocessors/AFMSublatticeSuperexchangeEnergy

## Overview

Computes the free energy in the simulation box (volume \Omega) due to the short-range superexchange interaction,

\begin{equation}
  \begin{aligned}
    F_\mathrm{super} = \int\limits_\Omega d^3 \mathbf{r} \left[ D_e \left(\mathbf{L}^2 - \mathbf{m}^2\right) \right]
  \end{aligned}
\end{equation}

where $D_e$ is the coupling coefficient driving collinearity of the two sublattices $\mathbf{m}_1$ and $\mathbf{m}_2$. The quantities $\mathbf{L}$ and $\mathbf{m}$ are the antiferromagnetic Néel and total magnetization vectors respectively. Their definitions are $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 - \mathbf{m}_2$.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMSublatticeSuperexchangeEnergy object.

!syntax parameters /Postprocessors/AFMSublatticeSuperexchangeEnergy

!syntax inputs /Postprocessors/AFMSublatticeSuperexchangeEnergy

!syntax children /Postprocessors/AFMSublatticeSuperexchangeEnergy
