# RotoBulkEnergyEighth

!syntax description /Postprocessors/RotoBulkEnergyEighth

## Overview

Computes the free energy in the simulation box (volume $\Omega$) corresponding to the bulk functional up to eighth order with,

\begin{equation}
  \begin{aligned}
    F_\mathrm{A,bulk} &= \int\limits_\Omega d^3 \mathbf{r} \left[ \beta__1 \left(A_x^2 + A_y^2 + A_z^2 \right) + \beta__{11} \left(A_x^4 + A_y^4 + A_z^4\right) + \beta__{12} \left(A_x^2 A_y^2 + A_y^2 A_z^2 + A_x^2 A_y^2\right)\right]. \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\beta__{111}\left(A_x^6 + A_y^6 + A_z^6\right) + \beta__{112}\left(A_x^2 \left(A_y^4 + A_z^4\right) + A_y^2 \left(A_x^4 + A_z^4\right)+A_z^2 \left(A_x^4 + A_y^4\right)\right)\right]. \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\beta__{123}\left(A_x^2 A_y^2 A_z^2\right) + \beta__{1111}\left(A_x^8 + A_y^8 + A_z^8\right) + \beta__{1122} \left(A_x^4 A_y^4 + A_y^4 A_z^4 + A_x^4 A_z^4\right)\right]. \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\beta__{1112}\left(A_x^6 \left(A_y^2 + A_z^2\right) + A_y^6 \left(A_x^2 + A_z^2\right) + A_z^6 \left(A_x^2 + A_y^2\right)\right)\right]. \\
    &+ \int\limits_\Omega d^3 \mathbf{r} \left[\beta__{1123}\left(A_x^4 A_y^2 A_z^2 + A_y^4 A_x^2 A_z^2 + A_z^4 A_x^2 A_y^2\right)\right]. \\
  \end{aligned}
\end{equation}

where $\mathbf{A}$ is the oxygen octahedral cage antiphase tilt vector and $\beta_{ij..}$ the Landau expansion coefficients of the phenomenological theory for a cubic parent phase.

## Example Input File Syntax

!! Describe and include an example of how to use the RotoBulkEnergyEighth object.

!syntax parameters /Postprocessors/RotoBulkEnergyEighth

!syntax inputs /Postprocessors/RotoBulkEnergyEighth

!syntax children /Postprocessors/RotoBulkEnergyEighth
