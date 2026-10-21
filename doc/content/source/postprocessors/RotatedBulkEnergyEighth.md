# RotatedBulkEnergyEighth

!syntax description /Postprocessors/RotatedBulkEnergyEighth

## Overview

Integrates the eighth-order Landau free energy of
[BulkEnergyEighth](BulkEnergyEighth.md) for a cubic crystal whose axes are turned about $z$ by
the material property `twist` (radians), i.e. the energy whose variation is
[RotatedBulkEnergyDerivativeEighth](RotatedBulkEnergyDerivativeEighth.md). The polynomial is
evaluated for the crystal-frame polarization $\mathbf{P}_c = R(\mathrm{twist})^T \mathbf{P}$.


## Example Input File Syntax

!listing test/tests/flexo/flexo_180DW_theta15.i block=Postprocessors/FbP

!syntax parameters /Postprocessors/RotatedBulkEnergyEighth

!syntax inputs /Postprocessors/RotatedBulkEnergyEighth

!syntax children /Postprocessors/RotatedBulkEnergyEighth
