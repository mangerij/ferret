# ComputeRotatedCubicParentElectrostrictiveStrain

!syntax description /Materials/ComputeRotatedCubicParentElectrostrictiveStrain

## Overview

Computes the electrostrictive eigenstrain $\varepsilon^Q_{ij} = Q_{ijkl} P_k P_l$ of a cubic parent
phase (with $\varepsilon^Q_{12} = Q_{44} P_1 P_2$, as in
[ComputeCubicParentElectrostrictiveStrain](ComputeCubicParentElectrostrictiveStrain.md)) whose crystal axes are turned
about $z$ by the material property `twist` (radians). At `twist` $= 0$ it reduces to the unrotated object.

!alert note title=Rotated crystal axes
The crystal axes are turned about $z$ by the material property `twist` (radians); the polarization,
its gradient and the strain are expressed in the crystal frame, the cubic expressions are evaluated
there and the result is turned back. For the displacement equation to see the same crystal,
`ComputeElasticityTensor` must be given `euler_angle_3` $= -$`twist` (in degrees).


## Example Input File Syntax

!listing test/tests/flexo/flexo_180DW_theta15.i block=Materials/electrostrictive_eigenstrain

!syntax parameters /Materials/ComputeRotatedCubicParentElectrostrictiveStrain

!syntax inputs /Materials/ComputeRotatedCubicParentElectrostrictiveStrain

!syntax children /Materials/ComputeRotatedCubicParentElectrostrictiveStrain
