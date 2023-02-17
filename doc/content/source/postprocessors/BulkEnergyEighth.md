# BulkEnergyEighth

!syntax description /Postprocessors/BulkEnergyEighth

## Overview

Computes the free energy in the simulation box (volume $\Omega$) corresponding to the bulk functional up to eighth order with,

\begin{equation}
  \begin{aligned}
    F_\mathrm{P,bulk} &= \int\limits_\Omega d^3 \mathbf{r} \left[ \alpha_1 \left(P_x^2 + P_y^2 + P_z^2 \right) + \alpha_{11} \left(P_x^4 + P_y^4 + P_z^4\right) + \alpha_{12} \left(P_x^2 P_y^2 + P_y^2 P_z^2 + P_x^2 P_y^2\right)\right] \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\alpha_{111}\left(P_x^6 + P_y^6 + P_z^6\right) + \alpha_{112}\left(P_x^2 \left(P_y^4 + P_z^4\right) + P_y^2 \left(P_x^4 + P_z^4\right)+P_z^2 \left(P_x^4 + P_y^4\right)\right)\right] \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\alpha_{123}\left(P_x^2 P_y^2 P_z^2\right) + \alpha_{1111}\left(P_x^8 + P_y^8 + P_z^8\right) + \alpha_{1122} \left(P_x^4 P_y^4 + P_y^4 P_z^4 + P_x^4 P_z^4\right)\right] \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\alpha_{1112}\left(P_x^6 \left(P_y^2 + P_z^2\right) + P_y^6 \left(P_x^2 + P_z^2\right) + P_z^6 \left(P_x^2 + P_y^2\right)\right)\right] \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\alpha_{1123}\left(P_x^4 P_y^2 P_z^2 + P_y^4 P_x^2 P_z^2 + P_z^4 P_x^2 P_y^2\right)\right] \\
  \end{aligned}
\end{equation}

where $\mathbf{P}$ is the polarization vector and $\alpha_{ij..}$ the Landau expansion coefficients of the phenomenological theory for a cubic parent phase.

## Example Input File Syntax

!! Describe and include an example of how to use the BulkEnergyEighth object.

!syntax parameters /Postprocessors/BulkEnergyEighth

!syntax inputs /Postprocessors/BulkEnergyEighth

!syntax children /Postprocessors/BulkEnergyEighth
