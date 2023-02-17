# BulkEnergyDensity

!syntax description /AuxKernels/BulkEnergyDensity

## Overview

Computes the free energy density corresponding to the bulk functional up to eighth order with,

\begin{equation}
  \begin{aligned}
    f_\mathrm{P,bulk} &= \alpha_1 \left(P_x^2 + P_y^2 + P_z^2 \right) + \alpha_{11} \left(P_x^4 + P_y^4 + P_z^4\right) + \alpha_{12} \left(P_x^2 P_y^2 + P_y^2 P_z^2 + P_x^2 P_y^2\right) \\
    &+ \alpha_{111}\left(P_x^6 + P_y^6 + P_z^6\right) + \alpha_{112}\left(P_x^2 \left(P_y^4 + P_z^4\right) + P_y^2 \left(P_x^4 + P_z^4\right)+P_z^2 \left(P_x^4 + P_y^4\right)\right) \\
    &+ \alpha_{123}\left(P_x^2 P_y^2 P_z^2\right) + \alpha_{1111}\left(P_x^8 + P_y^8 + P_z^8\right) + \alpha_{1122} \left(P_x^4 P_y^4 + P_y^4 P_z^4 + P_x^4 P_z^4\right) \\
    &+ \alpha_{1112}\left(P_x^6 \left(P_y^2 + P_z^2\right) + P_y^6 \left(P_x^2 + P_z^2\right) + P_z^6 \left(P_x^2 + P_y^2\right)\right) \\
    &+ \alpha_{1123}\left(P_x^4 P_y^2 P_z^2 + P_y^4 P_x^2 P_z^2 + P_z^4 P_x^2 P_y^2\right). \\
  \end{aligned}
\end{equation}

where $\mathbf{P}$ is the polarization vector and $\alpha_{ij..}$ the Landau expansion coefficients of the phenomenological theory for a cubic parent phase.


## Example Input File Syntax

!! Describe and include an example of how to use the BulkEnergyDensity object.

!syntax parameters /AuxKernels/BulkEnergyDensity

!syntax inputs /AuxKernels/BulkEnergyDensity

!syntax children /AuxKernels/BulkEnergyDensity
