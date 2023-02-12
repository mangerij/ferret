# PolarizationNWEMarker

!syntax description /Adaptivity/Markers/PolarizationNWEMarker

## Overview

Useful `Adaptivity` `Marker` allows for refinement and coarsening based on the values of the expected order parameter magnitude. Typically, in a ground state calculation (i.e. a film), the default values of the `Adaptivity` system will refine too much when $\mathbf{P}$ is near zero (the `ICs` of the calculation). This is a disadvantage because we are only concerned with refinment of the domain walls and coarsening of the homogeneous regions. The user can then control when the refinement and coarsening begins - which is dependent on the material's saturation polarization $P_s$.

## Example Input File Syntax

!! Describe and include an example of how to use the PolarizationNWEMarker object.

!syntax parameters /Adaptivity/Markers/PolarizationNWEMarker

!syntax inputs /Adaptivity/Markers/PolarizationNWEMarker

!syntax children /Adaptivity/Markers/PolarizationNWEMarker
