# FluctuationKernel

!alert construction title=Undocumented Class
The FluctuationKernel has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Kernels/FluctuationKernel

## Overview

Introduces a residual contribution of the form

\begin{equation}
  \begin{aligned}
    -\left(\psi_h, \Delta \pi \right) = 0
  \end{aligned}
\end{equation}

where $\psi_h$ is the test function and $\Delta \pi$ is a small number. This can be applied to any variable to introduce random fluctuations on the order of $\Delta \pi$ which can be useful in quasi-static hysteresis calculations.

## Example Input File Syntax

!! Describe and include an example of how to use the FluctuationKernel object.

!syntax parameters /Kernels/FluctuationKernel

!syntax inputs /Kernels/FluctuationKernel

!syntax children /Kernels/FluctuationKernel
