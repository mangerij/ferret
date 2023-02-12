# ScalarBulkEnergyA

!syntax description /ScalarKernels/ScalarBulkEnergyA

## Overview

Sets up the RHS of the coupled Landau-Khalatnikov equation. Since we are working with `ScalarKernels`, this is effectively a single grid point, or homogeneous calculation. In the context of micromagnetics, this is the macrospin analog for ferroelectrics. This object concerns the bulk contribution from the antiphase oxygen octahedral tilts $\mathbf{A}$ (in the case of $\mathrm{BiFeO}_3$),

\begin{equation}
  \begin{aligned}
     \frac{\partial \mathbf{A}}{\partial t} = - \Gamma_A \frac{\delta f_\mathrm{A,bulk}}{\delta \mathbf{A}}
  \end{aligned}
\end{equation}

Ignoring the LHS and working in index notation,

\begin{equation}
  \begin{aligned}
     - \Gamma_A \frac{\delta f_\mathrm{A,bulk}}{\delta A_\eta} &= - \Gamma_A \frac{\partial}{\partial A_\eta} \left( \beta_{ij} A_i A_j + \beta_{ijkl} A_i A_j A_k A_l\right) \\
     &= - \Gamma_A \left( 2 \beta_{i\eta} A_i  + 4 \beta_{ijk\eta} A_i A_j A_k\right) \\
  \end{aligned}
\end{equation}

In Voight notation, many of these terms vanish due to symmetry of the parent phase.

## Example Input File Syntax

!! Describe and include an example of how to use the ScalarBulkEnergyA object.

!syntax parameters /ScalarKernels/ScalarBulkEnergyA

!syntax inputs /ScalarKernels/ScalarBulkEnergyA

!syntax children /ScalarKernels/ScalarBulkEnergyA
