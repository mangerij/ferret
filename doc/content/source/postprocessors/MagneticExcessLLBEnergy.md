# MagneticExcessLLBEnergy

!syntax description /Postprocessors/MagneticExcessLLBEnergy

## Overview

Useful tracker of the Landau-Lifshitz-Bloch implementation where,

\begin{equation}
  \begin{aligned}
    F_\mathrm{excess} = \int\limits_\Omega d^3 \mathbf{r} \left\{\frac{1}{4} \left(\frac{\tilde{\alpha}_{LLB}}{1+\alpha^2}\right)\left(|\mathbf{m}| - 1\right) \right\}
  \end{aligned}
\end{equation}

is computed in the simulation box (volume $\Omega$). This `Postprocessor` tracks the deviations of the normalized magnetization vector $\mathbf{m}$ from unity. The energy here is arbitrary (and unitless) but can be cast into real units by using the prefactor `energy`$\_$`scale` included in the object input parameters.

## Example Input File Syntax

!! Describe and include an example of how to use the MagneticExcessLLBEnergy object.

!syntax parameters /Postprocessors/MagneticExcessLLBEnergy

!syntax inputs /Postprocessors/MagneticExcessLLBEnergy

!syntax children /Postprocessors/MagneticExcessLLBEnergy
