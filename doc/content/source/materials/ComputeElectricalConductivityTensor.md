# ComputeElectricalConductivityTensor

!syntax description /Materials/ComputeElectricalConductivityTensor

## Overview

Computes (and stores) the electrical conductivity tensor $\sigma_{ij}$ at each quadrature point in the finite element mesh. The orientation of the tensor can be rotated via,

\begin{equation}
  \begin{aligned}
    \tilde{\sigma}_{ij} = R_{i\beta}R_{j\gamma} \sigma_{\beta\gamma}
  \end{aligned}
\end{equation}

for some arbitrary rotation matrices $R_{ij}$ using the internal `RotationTensor` operation in MOOSE utils. The rotation operator $R_{ij}$ accepts Euler angles in the standard Bunge sequence ($\mathbf{ZXZ}$). Note that $\sigma_{ij}$ is generally the nomenclature for electrical conductivty, not to be confused with the linear elastic stress.

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeElectricalConductivityTensor object.

!syntax parameters /Materials/ComputeElectricalConductivityTensor

!syntax inputs /Materials/ComputeElectricalConductivityTensor

!syntax children /Materials/ComputeElectricalConductivityTensor
