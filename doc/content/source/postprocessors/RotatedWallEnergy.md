# RotatedWallEnergy

!syntax description /Postprocessors/RotatedWallEnergy

## Overview

Integrates the gradient energy $f_\mathrm{grad} = \frac{1}{2} G_{ijkl}\, \partial_j P_i\, \partial_l P_k$
of [RotatedWallEnergyDerivative](RotatedWallEnergyDerivative.md) (bare coefficients `G11`, `G12`,
`G44`) for a cubic crystal whose axes are turned about $z$ by the material property `twist`.


## Example Input File Syntax

!listing test/tests/flexo/flexo_180DW_theta15.i block=Postprocessors/FgP

!syntax parameters /Postprocessors/RotatedWallEnergy

!syntax inputs /Postprocessors/RotatedWallEnergy

!syntax children /Postprocessors/RotatedWallEnergy
