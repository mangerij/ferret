# ComputeGCoeffTensor

!syntax description /Materials/ComputeGCoeffTensor

## Overview

Computes (and stores) a polar-optic coefficient tensor at each quadrature point of the finite element mesh which can be rotated,

\begin{equation}
  \begin{aligned}
    \tilde{g}_{ijkl} = R_{i\beta} R_{j\gamma} R_{k\delta} R_{l\eta} g_{\beta\gamma\delta\eta}
  \end{aligned}
\end{equation}

with some arbitrary rotation matrices $R_{ij}$.


## Example Input File Syntax

!! Describe and include an example of how to use the ComputeGCoeffTensor object.

!syntax parameters /Materials/ComputeGCoeffTensor

!syntax inputs /Materials/ComputeGCoeffTensor

!syntax children /Materials/ComputeGCoeffTensor
