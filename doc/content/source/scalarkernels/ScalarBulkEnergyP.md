# ScalarBulkEnergyP

!syntax description /ScalarKernels/ScalarBulkEnergyP

## Overview

Sets up the RHS of the coupled Landau-Khalatnikov equation. Since we are working with `ScalarKernels`, this is effectively a single grid point, or homogeneous calculation. In the context of micromagnetics, this is the macrospin analog for ferroelectrics. This object concerns the bulk contribution from the electric polarization $\mathbf{P}$ (in the case of $\mathrm{BiFeO}_3$),

\begin{equation}
  \begin{aligned}
     \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta f_\mathrm{P,bulk}}{\delta \mathbf{P}}
  \end{aligned}
\end{equation}

Ignoring the LHS and working in index notation,

\begin{equation}
  \begin{aligned}
     - \Gamma_P \frac{\delta f_\mathrm{P,bulk}}{\delta P_\eta} &= - \Gamma_P \frac{\partial}{\partial P_\eta} \left( \alpha_{ij} P_i P_j + \alpha_{ijkl} P_i P_j P_k P_l\right) \\
     &= - \Gamma_P \left( 2 \alpha_{i\eta} P_i  + 4 \alpha_{ijk\eta} P_i P_j P_k\right) \\
  \end{aligned}
\end{equation}

In Voight notation, many of these terms vanish due to symmetry of the parent phase.

## Example Input File Syntax

!! Describe and include an example of how to use the ScalarBulkEnergyP object.

!syntax parameters /ScalarKernels/ScalarBulkEnergyP

!syntax inputs /ScalarKernels/ScalarBulkEnergyP

!syntax children /ScalarKernels/ScalarBulkEnergyP
