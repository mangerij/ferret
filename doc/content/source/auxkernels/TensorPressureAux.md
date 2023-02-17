# TensorPressureAux

!syntax description /AuxKernels/TensorPressureAux

## Overview

Calculates the trace of the stress tensor $\sigma_{ij}$ (hydrostatic pressure) according to

\begin{equation}
  \begin{aligned}
    p &= - \frac{1}{3} tr(\sigma_{ij}) \\
      &= - \frac{1}{3} \left\{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}\right\} \\
  \end{aligned}
\end{equation}

## Example Input File Syntax

!! Describe and include an example of how to use the TensorPressureAux object.

!syntax parameters /AuxKernels/TensorPressureAux

!syntax inputs /AuxKernels/TensorPressureAux

!syntax children /AuxKernels/TensorPressureAux
