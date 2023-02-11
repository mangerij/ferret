# WallEnergyDensity

!syntax description /AuxKernels/WallEnergyDensity

## Overview

Computes the free energy density due to gradients in the polarization vector field $\mathbf{P}$. It is defined as follows,

\begin{equation}
  \begin{aligned}&f_\mathrm{\nabla P} =\frac{G_{11}}{2}  \left( P_{x,x}^2 + P_{y,y}^2 + P_{z,z}^2 \right) +  G_{12}  \left(P_{x,x} P_{y,y} + P_{y,y} P_{z,z} + P_{x,x} P_{z,z} \right) + \frac{G_{44}}{2} \left[\left(P_{x,y} + P_{y,x} \right)^2+ \left(P_{y,z} + P_{z,y} \right)^2 + \left(P_{x,z} + P_{z,x}\right)^2\right]
  \end{aligned}
\end{equation}

with gradient coefficients $G_{ij}$. Here $P_{i,j}$ denotes the $i^\mathrm{th}$ component and $j^\mathrm{th}$ derivative direction.

## Example Input File Syntax

!! Describe and include an example of how to use the WallEnergyDensity object.

!syntax parameters /AuxKernels/WallEnergyDensity

!syntax inputs /AuxKernels/WallEnergyDensity

!syntax children /AuxKernels/WallEnergyDensity
