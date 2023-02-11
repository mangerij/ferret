# AFMExchangeStiffnessEnergy

!syntax description /Postprocessors/AFMExchangeStiffnessEnergy

## Overview

Computes the free energy in the simulation box (volume $\Omega$) due to the long-range exchange interactions in the antiferromagnet,

\begin{equation}
  \begin{aligned}
    F_\mathrm{\nabla L} = \int\limits_\Omega \left\{ A_e \left[\left(\nabla L_x\right)^2 + \left(\nabla L_y\right)^2 \left(\nabla L_z\right)^2\right] \right\}.
  \end{aligned}
\end{equation}

The quantity $\mathbf{L}$ is the antiferromagnetic Néel vector along with coupling constant $A_e$.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMExchangeStiffnessEnergy object.

!syntax parameters /Postprocessors/AFMExchangeStiffnessEnergy

!syntax inputs /Postprocessors/AFMExchangeStiffnessEnergy

!syntax children /Postprocessors/AFMExchangeStiffnessEnergy
