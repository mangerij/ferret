# BandGapAuxZnO

!syntax description /AuxKernels/BandGapAuxZnO

## Overview

Calculates a simple band gap shift due to stress tensor components $\sigma_{xx}$ and $\sigma_{zz}$,

\begin{equation}
  \begin{aligned}
    E_g = E_g^0 + \frac{1}{2}\frac{1}{1-R_b}\left\[\left(d_b + d_u R_b\right) \varepsilon_{xx} + \left(d_u + \nu d_b\right)\varepsilon_{zz}\right\]
  \end{aligned}
\end{equation}

where $E_g^0$ is the unstrained band gap energy and the other parameters correspond to the mechanical deformation.

## Example Input File Syntax

!! Describe and include an example of how to use the BandGapAuxZnO object.

!syntax parameters /AuxKernels/BandGapAuxZnO

!syntax inputs /AuxKernels/BandGapAuxZnO

!syntax children /AuxKernels/BandGapAuxZnO
