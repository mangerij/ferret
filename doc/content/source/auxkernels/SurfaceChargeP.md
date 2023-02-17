# SurfaceChargeP

!syntax description /AuxKernels/SurfaceChargeP

## Overview

Computes the surface charge due to the polarization $\mathbf{P}$ which is typically denoted by $\sigma_b$

\begin{equation}
  \begin{aligned}
    \sigma_b = \mathbf{P}\cdot\mathbf{n}
  \end{aligned}
\end{equation}

where $\mathbf{n}$ is a surface normal. This expression can be computed over arbitrary curved surfaces.

## Example Input File Syntax

!! Describe and include an example of how to use the SurfaceChargeP object.

!syntax parameters /AuxKernels/SurfaceChargeP

!syntax inputs /AuxKernels/SurfaceChargeP

!syntax children /AuxKernels/SurfaceChargeP
