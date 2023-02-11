# ComputeElectrostrictiveTensor

!alert construction title=Undocumented Class
The ComputeElectrostrictiveTensor has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Materials/ComputeElectrostrictiveTensor

## Overview

Calculates the electrostrictive tensor components $q_{ijkl}$. Here,

\begin{equation}
  \begin{aligned}
      q_{ijkl} = C_{ijmn} Q_{mnkl}.
  \end{aligned}
\end{equation}

with $C_{ijmn}$ and $Q_{mnkl}$ are the elastic stiffness and electrostictive coefficients respectively. The $q_{ijkl}$ tensor may be rotated via,

\begin{equation}
  \begin{aligned}
    \tilde{q}_{ijkl} = R_{i\alpha} R_{j\beta} R_{k \gamma} R_{l\delta} q_{\alpha\beta\gamma\delta}.
  \end{aligned}
\end{equation}

with the internal `RotationTensor` operation in MOOSE utils. The rotation operator $R_{ij}$ accepts Euler angles in the standard Bunge sequence ($\mathbf{ZXZ}$).

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeElectrostrictiveTensor object.

!syntax parameters /Materials/ComputeElectrostrictiveTensor

!syntax inputs /Materials/ComputeElectrostrictiveTensor

!syntax children /Materials/ComputeElectrostrictiveTensor
