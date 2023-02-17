# ChangeInRefractiveIndexWithGCoeffPolar

!syntax description /AuxKernels/ChangeInRefractiveIndexWithGCoeffPolar

## Overview

Computes the change in the refractive indices,

\begin{equation}
  \begin{aligned}
    \Delta n_{ij} = \frac{1}{2} B_{ij}^3 \left[\Delta B_{kl} + \Delta B_{kl}(\mathbf{P})\right].
  \end{aligned}
\end{equation}

where $\Delta B_{ij}$ is the change in the indicatrix (i.e. due to elastooptic effects) and $\Delta B_{ij}(\mathbf{P})$ is due to a polar phase transition. The above expression does not use index summation so typically one picks $i = k$ and $j = l$. NOTE: unsure why we have a minus sign here, will change this soon.

## Example Input File Syntax

!! Describe and include an example of how to use the ChangeInRefractiveIndexWithGCoeffPolar object.

!syntax parameters /AuxKernels/ChangeInRefractiveIndexWithGCoeffPolar

!syntax inputs /AuxKernels/ChangeInRefractiveIndexWithGCoeffPolar

!syntax children /AuxKernels/ChangeInRefractiveIndexWithGCoeffPolar
