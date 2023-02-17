# ChangeInRefractiveIndexElectro

!syntax description /AuxKernels/ChangeInRefractiveIndexElectro

## Overview

Computes the change in the refractive indices,

\begin{equation}
  \begin{aligned}
    \Delta n_{ij} = - \frac{1}{2} B_{ij}^3 \Delta B_{kl}(\mathbf{E}).
  \end{aligned}
\end{equation}

where $\Delta B_{ij}(\mathbf{E})$ is the change in the indicatrix due an electric field $\mathbf{E}$. The above expression does not use index summation so typically one picks $i = k$ and $j = l$.


## Example Input File Syntax

!! Describe and include an example of how to use the ChangeInRefractiveIndexElectro object.

!syntax parameters /AuxKernels/ChangeInRefractiveIndexElectro

!syntax inputs /AuxKernels/ChangeInRefractiveIndexElectro

!syntax children /AuxKernels/ChangeInRefractiveIndexElectro
