# ComputeSeebeckTDepTensor

!syntax description /Materials/ComputeSeebeckTDepTensor

## Overview

Computes a temperature dependent Seebeck tensor $S_{ij}(T)$,

\begin{equation}
  \begin{aligned}
    S_{ij}(T) = a_{ij}^S + b_{ij}^S T + c_{ij}^S T^2
  \end{aligned}
\end{equation}

where $T$ is the temperature and $a_{ij}^S, b_{ij}^S$, and $c_{ij}^S$ are coupling tensor coefficients of rank four. They are seeded via `fillFromInputVector` which requires nine entries.

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeSeebeckTDepTensor object.

!syntax parameters /Materials/ComputeSeebeckTDepTensor

!syntax inputs /Materials/ComputeSeebeckTDepTensor

!syntax children /Materials/ComputeSeebeckTDepTensor
