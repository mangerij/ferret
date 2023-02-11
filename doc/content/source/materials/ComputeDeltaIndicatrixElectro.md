# ComputeDeltaIndicatrixElectro

!alert construction title=Undocumented Class
The ComputeDeltaIndicatrixElectro has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Materials/ComputeDeltaIndicatrixElectro

## Overview

Computes the change of the indicatrix components $\Delta B_{ij}$ due to,

\begin{equation}
  \begin{aligned}
    \Delta B_{ij} = r_{ijk} E_k,
  \end{aligned}
\end{equation}

where $r_{ijk}$ is the electrooptic tensor of rank three. The variable $E_k$ is the component of the electric field related to the electrostatic potential $\Phi_\mathrm{E}$ via $\mathbf{E} = - \nabla \Phi_\mathrm{E}$ which is calculated with the appropriate Laplace or Poisson equation `Kernels`.

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeDeltaIndicatrixElectro object.

!syntax parameters /Materials/ComputeDeltaIndicatrixElectro

!syntax inputs /Materials/ComputeDeltaIndicatrixElectro

!syntax children /Materials/ComputeDeltaIndicatrixElectro
