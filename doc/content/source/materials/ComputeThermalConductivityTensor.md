# ComputeThermalConductivityTensor

!syntax description /Materials/ComputeThermalConductivityTensor

## Overview

Computes (and stores) the thermal conductivity tensor $\kappa_{ij}$ at each quadrature point in the finite element mesh. The orientation of the tensor can be rotated via,

\begin{equation}
  \begin{aligned}
    \tilde{\kappa}_{ij} = R_{i\beta}R_{j\gamma} \kappa_{\beta\gamma}
  \end{aligned}
\end{equation}

for some arbitrary rotation matrices $R_{ij}$ using the internal `RotationTensor` operation in MOOSE utils. The rotation operator $R_{ij}$ accepts Euler angles in the standard Bunge sequence ($\mathbf{ZXZ}$).

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeThermalConductivityTensor object.

!syntax parameters /Materials/ComputeThermalConductivityTensor

!syntax inputs /Materials/ComputeThermalConductivityTensor

!syntax children /Materials/ComputeThermalConductivityTensor
