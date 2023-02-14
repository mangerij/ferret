# RotostrictiveCouplingDispDerivative

!syntax description /Kernels/RotostrictiveCouplingDispDerivative

## Overview

We now consider the condition for mechanical equilibrium. In addition to the electric polarization, $\mathbf{P}$, there are ferroelectrics where secondary structural order parameters condense below the phase transition temperature. In the case of $\mathrm{BiFeO}_3$, the oxygen octahedral tilt cages can generate a contribution to the total spontaneous strain,

\begin{equation}
  \begin{aligned}
    \varepsilon_{ij}^\mathrm{eig,A} = R_{ijkl} A_k A_l,
  \end{aligned}
\end{equation}

where $R_{ijkl}$ is the rotostrictive coefficient tensor and $\mathbf{A}$ is the order parameter associated with the antiphase tilts.

## Example Input File Syntax

!! Describe and include an example of how to use the RotostrictiveCouplingDispDerivative object.

!syntax parameters /Kernels/RotostrictiveCouplingDispDerivative

!syntax inputs /Kernels/RotostrictiveCouplingDispDerivative

!syntax children /Kernels/RotostrictiveCouplingDispDerivative
