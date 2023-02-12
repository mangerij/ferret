# ComputeThermalConductivityTDepTensor

!syntax description /Materials/ComputeThermalConductivityTDepTensor

## Overview

Computes a temperature dependent thermal conductivity tensor $\kappa_{ij}(T)$,

\begin{equation}
  \begin{aligned}
    \kappa_{ij}(T) = a_{ij}^\kappa + b_{ij}^\kappa T + c_{ij}^\kappa T^2
  \end{aligned}
\end{equation}

where $T$ is the temperature and $a_{ij}^\kappa, b_{ij}^\kappa$, and $c_{ij}^\kappa$ are coupling tensor coefficients of rank four. They are seeded via `fillFromInputVector` which requires nine entries.

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeThermalConductivityTDepTensor object.

!syntax parameters /Materials/ComputeThermalConductivityTDepTensor

!syntax inputs /Materials/ComputeThermalConductivityTDepTensor

!syntax children /Materials/ComputeThermalConductivityTDepTensor
