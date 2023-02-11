# ComputeElectroopticTensor

!alert construction title=Undocumented Class
The ComputeElectroopticTensor has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Materials/ComputeElectroopticTensor

## Overview

Stores the electrooptic tensor components $r_{ijk}$. The tensor may be rotated via,

\begin{equation}
  \begin{aligned}
    \tilde{r}_{ijk} = R_{i\alpha} R_{j\beta} R_{k\gamma} r_{\alpha\beta\gamma}.
  \end{aligned}
\end{equation}

with the internal `RotationTensor` operation in MOOSE utils. The rotation operator $R_{ij}$ accepts Euler angles in the standard Bunge sequence ($\mathbf{ZXZ}$).


## Example Input File Syntax

!! Describe and include an example of how to use the ComputeElectroopticTensor object.

!syntax parameters /Materials/ComputeElectroopticTensor

!syntax inputs /Materials/ComputeElectroopticTensor

!syntax children /Materials/ComputeElectroopticTensor
