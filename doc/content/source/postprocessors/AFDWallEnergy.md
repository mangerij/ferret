# AFDWallEnergy

!syntax description /Postprocessors/AFDWallEnergy

## Overview

Computes the free energy in the simulation box (volume $\Omega$) due to gradients in the antiphase oxygen octahedral tilt field $\mathbf{A}$. It is defined as follows,

\begin{equation}
  \begin{aligned}
  F_\mathrm{\nabla A} = \int\limits_\Omega d^3 \mathbf{r} \left\{ \frac{H_{11}}{2}   \left( A_{x,x}^2 + A_{y,y}^2 + A_{z,z}^2 \right) +  H_{12}  \left(A_{x,x} A_{y,y} + A_{y,y} A_{z,z} + A_{x,x} A_{z,z} \right) + \frac{H_{44}}{2} \left[\left(A_{x,y} + A_{y,x} \right)^2+ \left(A_{y,z} + A_{z,y} \right)^2 + \left(A_{x,z} + A_{z,x}\right)^2\right] \right\}
  \end{aligned}
\end{equation}

with gradient coefficients $H_{ij}$. Here $A_{i,j}$ denotes the $i^\mathrm{th}$ component and $j^\mathrm{th}$ derivative direction.

## Example Input File Syntax

!! Describe and include an example of how to use the AFDWallEnergy object.

!syntax parameters /Postprocessors/AFDWallEnergy

!syntax inputs /Postprocessors/AFDWallEnergy

!syntax children /Postprocessors/AFDWallEnergy
