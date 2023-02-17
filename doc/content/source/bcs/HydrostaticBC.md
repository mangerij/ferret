# HydrostaticBC

!syntax description /BCs/HydrostaticBC

## Overview

Computes a residual contribution using the `IntegratedBC` class with $p$ a specified hydrostatic scalar pressure,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{u_i} = -\left\langle\psi_h, p \frac{u_i}{|\mathbf{u}|} \cdot \mathbf{n}\right\rangle
  \end{aligned}
\end{equation}

where $\mathbf{u}$ is the elastic displacement vector, $\langle, \rangle$ is a surface integral, $\psi_h$ a test function of the finite element method, $\mathbf{n}$ the surface normal.

## Example Input File Syntax

!! Describe and include an example of how to use the HydrostaticBC object.

!syntax parameters /BCs/HydrostaticBC

!syntax inputs /BCs/HydrostaticBC

!syntax children /BCs/HydrostaticBC
