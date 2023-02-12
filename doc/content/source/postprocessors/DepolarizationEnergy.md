# DepolarizationEnergy

!syntax description /Postprocessors/DepolarizationEnergy

## Overview

Computes a volume integral over the free energy density defined by,

\begin{equation}
  \begin{aligned}
     F_\mathrm{depol} = \frac{1}{2} \int\limits_\Omega d^3 \mathbf{r} \lambda \frac{1}{\epsilon_b} P_z \bar{P}_z
  \end{aligned}
\end{equation}

with $\lambda$ a screening length. The value $\bar{P}_z$ is a volume averaged quantity computed from `ElementAverageValue`. This piece of code is restrictive in the sense that it is hardcoded for only the $\hat{z}$ direction. It can be generalized further if needed.

## Example Input File Syntax

!! Describe and include an example of how to use the DepolarizationEnergy object.

!syntax parameters /Postprocessors/DepolarizationEnergy

!syntax inputs /Postprocessors/DepolarizationEnergy

!syntax children /Postprocessors/DepolarizationEnergy
