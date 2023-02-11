# ChangeInRefractiveIndex

!alert construction title=Undocumented Class
The ChangeInRefractiveIndex has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/ChangeInRefractiveIndex

## Overview

Computes the change in the refractive indices,

\begin{equation}
  \begin{aligned}
    \Delta n_{ij} = - \frac{1}{2} B_{ij}^3 \Delta B_{kl}.
  \end{aligned}
\end{equation}

where $\Delta B_{ij}$ is the change in the indicatrix. The above expression does not use index summation so typically one picks $i = k$ and $j = l$.


## Example Input File Syntax

!! Describe and include an example of how to use the ChangeInRefractiveIndex object.

!syntax parameters /AuxKernels/ChangeInRefractiveIndex

!syntax inputs /AuxKernels/ChangeInRefractiveIndex

!syntax children /AuxKernels/ChangeInRefractiveIndex
