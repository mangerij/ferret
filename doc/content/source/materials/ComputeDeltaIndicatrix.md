# ComputeDeltaIndicatrix

!alert construction title=Undocumented Class
The ComputeDeltaIndicatrix has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Materials/ComputeDeltaIndicatrix

## Overview

Computes the change of the indicatrix components $\Delta B_{ij}$ due to,

\begin{equation}
  \begin{aligned}
    \Delta B_{ij} = p_{ijkl} \varepsilon_{kl},
  \end{aligned}
\end{equation}

where $p_{ijkl}$ is the elastooptic tensor related to the photoelastic tensor via $\pi_{ijmn}$ via the relationship

\begin{equation}
  \begin{aligned}
    p_{ijkl} = \pi_{ijmn} C_{mnkl}
  \end{aligned}
\end{equation}

with $C_{mnkl}$ the elastic stiffness tensor. The variable $\varepsilon_{ij}$ is the linear elastic strain tensor as computed by `TensorMechanics`.

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeDeltaIndicatrix object.

!syntax parameters /Materials/ComputeDeltaIndicatrix

!syntax inputs /Materials/ComputeDeltaIndicatrix

!syntax children /Materials/ComputeDeltaIndicatrix
