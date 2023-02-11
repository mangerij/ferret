# ElectronCurrentDensityBC

!alert construction title=Undocumented Class
The ElectronCurrentDensityBC has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /BCs/ElectronCurrentDensityBC

## Overview

Sets an `IntegratedBC` condition on a surface providing a residual term of the form,

\begin{equation}
  \begin{aligned}
    R_u = -q \mu_n N_c e^{\frac{q u - E_c}{k_B T}} \nabla u \cdot \mathbf{n} - \frac{k_B T}{q} N_c e^{\frac{q \nabla u \cdot \mathbf{n} - E_c}{k_B T} \mathbf{n}.
  \end{aligned}
\end{equation}

for some variable $u$ (usually the electrostatic potential).

## Example Input File Syntax

!! Describe and include an example of how to use the ElectronCurrentDensityBC object.

!syntax parameters /BCs/ElectronCurrentDensityBC

!syntax inputs /BCs/ElectronCurrentDensityBC

!syntax children /BCs/ElectronCurrentDensityBC
