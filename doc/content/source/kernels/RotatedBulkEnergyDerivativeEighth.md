# RotatedBulkEnergyDerivativeEighth

!syntax description /Kernels/RotatedBulkEnergyDerivativeEighth

## Overview

Calculates the residual contribution of the eighth-order Landau free energy of
[BulkEnergyDerivativeEighth](BulkEnergyDerivativeEighth.md) for a cubic crystal whose axes are
turned about $z$ by the material property `twist` (radians): the Landau polynomial is evaluated
for the crystal-frame polarization $\mathbf{P}_c = R(\mathrm{twist})^T \mathbf{P}$. The energy is
computed by [RotatedBulkEnergyEighth](RotatedBulkEnergyEighth.md).

!alert note title=Rotated crystal axes
The crystal axes are turned about $z$ by the material property `twist` (radians); the polarization,
its gradient and the strain are expressed in the crystal frame, the cubic expressions are evaluated
there and the result is turned back. For the displacement equation to see the same crystal,
`ComputeElasticityTensor` must be given `euler_angle_3` $= -$`twist` (in degrees).


## Example Input File Syntax

!listing test/tests/flexo/flexo_180DW_theta15.i block=Kernels/bed_x

!syntax parameters /Kernels/RotatedBulkEnergyDerivativeEighth

!syntax inputs /Kernels/RotatedBulkEnergyDerivativeEighth

!syntax children /Kernels/RotatedBulkEnergyDerivativeEighth
