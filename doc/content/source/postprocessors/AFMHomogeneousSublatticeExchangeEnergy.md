# AFMHomogeneousSublatticeExchangeEnergy

!syntax description /Postprocessors/AFMHomogeneousSublatticeExchangeEnergy

## Overview

Computes the free energy in the simulation box (volume $\Omega$) due to the short-range sublattice exchange energy between two sublattice magnetizations $\{\mathbf{m}_1, \mathbf{m}_2\}$ in an antiferromagnet corresponding to,

\begin{equation}
  \begin{aligned}
    F_\mathrm{exch} &= \int\limits_\Omega d^3 \mathbf{r} \, 4 D_e \left[ \mathbf{L}^2 - \mathbf{m}^2 \right],
  \end{aligned}
\end{equation}

where $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$. Note that we do not use the conventional factor of 2 in these definitions but they can be evaluated after or during the simulation.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMHomogeneousSublatticeExchangeEnergy object.

!syntax parameters /Postprocessors/AFMHomogeneousSublatticeExchangeEnergy

!syntax inputs /Postprocessors/AFMHomogeneousSublatticeExchangeEnergy

!syntax children /Postprocessors/AFMHomogeneousSublatticeExchangeEnergy
