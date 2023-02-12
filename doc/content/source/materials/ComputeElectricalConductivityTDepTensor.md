# ComputeElectricalConductivityTDepTensor

!syntax description /Materials/ComputeElectricalConductivityTDepTensor

## Overview

Computes a temperature dependent electrical conductivity tensor $\sigma_{ij}(T)$,

\begin{equation}
  \begin{aligned}
    \sigma_{ij}(T) = a_{ij}^\sigma + b_{ij}^\sigma T + c_{ij}^\sigma T^2
  \end{aligned}
\end{equation}

where $T$ is the temperature and $a_{ij}^\sigma, b_{ij}^\sigma$, and $c_{ij}^\sigma$ are coupling tensor coefficients of rank four. They are seeded via `fillFromInputVector` which requires nine entries. Note that $\sigma_{ij}$ is generally the nomenclature for electrical conductivty, not to be confused with the linear elastic stress.

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeElectricalConductivityTDepTensor object.

!syntax parameters /Materials/ComputeElectricalConductivityTDepTensor

!syntax inputs /Materials/ComputeElectricalConductivityTDepTensor

!syntax children /Materials/ComputeElectricalConductivityTDepTensor
