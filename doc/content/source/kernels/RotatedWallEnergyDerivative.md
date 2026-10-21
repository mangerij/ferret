# RotatedWallEnergyDerivative

!syntax description /Kernels/RotatedWallEnergyDerivative

## Overview

Calculates the residual contribution of the gradient (domain wall) energy

\begin{equation}
  f_\mathrm{grad} = \frac{1}{2} G_{ijkl} \frac{\partial P_i}{\partial x_j}\frac{\partial P_k}{\partial x_l}
\end{equation}

with the cubic tensor $G_{1111} = G_{11}$, $G_{1122} = G_{12}$, $G_{1212} = G_{44}$, for a crystal
whose axes are turned about $z$ by the material property `twist` (radians). Unlike
[WallEnergyDerivative](WallEnergyDerivative.md), which takes `G110` and the ratios `G11_G110`, ...,
this object takes the bare coefficients `G11`, `G12`, `G44` and has no antisymmetric (`G44P`) term;
at `twist` $= 0$ the two agree for `G110 = G11`, `G11_G110 = 1`, `G12_G110 = G12/G11`,
`G44_G110 = G44/G11`, `G44P_G110 = 0`. The energy is computed by [RotatedWallEnergy](RotatedWallEnergy.md).


## Example Input File Syntax

!listing test/tests/flexo/flexo_180DW_theta15.i block=Kernels/walled_x

!syntax parameters /Kernels/RotatedWallEnergyDerivative

!syntax inputs /Kernels/RotatedWallEnergyDerivative

!syntax children /Kernels/RotatedWallEnergyDerivative
