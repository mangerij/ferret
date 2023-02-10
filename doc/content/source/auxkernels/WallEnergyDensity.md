# WallEnergyDensity

!alert construction title=Undocumented Class
The WallEnergyDensity has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/WallEnergyDensity

## Overview

Computes the local free energy density due to gradients in $\mathbf{P}$,

\begin{equation}
  \begin{aligned}
  &f_\mathrm{\nabla P} = \frac{G_{11}}{2}  \left( P_{x,x}^2 + P_{y,y}^2 + P_{z,z}^2 \right) \\
  &+  G_{12}  \left(P_{x,x} P_{y,y} + P_{y,y} P_{z,z} + P_{x,x} P_{z,z} \right) \\ 
  &+ \frac{G_{44}}{2} \left[\left(P_{x,y} + P_{y,x} \right)^2+ \left(P_{y,z} + P_{z,y} \right)^2 + \left(P_{x,z} + P_{z,x}\right)^2\right].
  \end{aligned}
\end{equation}

## Example Input File Syntax

!! Describe and include an example of how to use the WallEnergyDensity object.

!syntax parameters /AuxKernels/WallEnergyDensity

!syntax inputs /AuxKernels/WallEnergyDensity

!syntax children /AuxKernels/WallEnergyDensity
