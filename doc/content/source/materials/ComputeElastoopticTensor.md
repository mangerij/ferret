# ComputeElastoopticTensor

!alert construction title=Undocumented Class
The ComputeElastoopticTensor has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Materials/ComputeElastoopticTensor

## Overview

Stores the elastooptic tensor components $p_{ijkl}$. The tensor may be rotated via,

\begin{equation}
  \begin{aligned}
    \tilde{p}_{ijkl} = R_{i\alpha} R_{j\beta} R_{k \gamma} R_{l\delta} p_{\alpha\beta\gamma\delta}.
  \end{aligned}
\end{equation}

with the internal `RotationTensor` operation in MOOSE utils. The rotation operator $R_{ij}$ accepts Euler angles in the standard Bunge sequence ($\mathbf{ZXZ}$).

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeElastoopticTensor object.

!syntax parameters /Materials/ComputeElastoopticTensor

!syntax inputs /Materials/ComputeElastoopticTensor

!syntax children /Materials/ComputeElastoopticTensor
