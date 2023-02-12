# ThermoelectricMaterial

!syntax description /Materials/ThermoelectricMaterial

## Overview

Uses the automatic differentation (AD) system to compute derivatives with respect to temperature $T$ of the materials coefficients $S_{ij}(T)$, $\sigma_{ij}(T)$, and $\kappa_{ij}(T)$. These tensors are defined in `ComputeSeebeckTDepTensor`, `ComputeElectricalConductivityTDepTensor`, and `ComputeThermalConductivityTDepTensor`. We should emphasize that typically in the literature, the nomenclature for the electrical conductivity tensor is $\sigma_{ij}$ but this is not to be confused with the linear elastic stress tensor.

## Example Input File Syntax

!! Describe and include an example of how to use the ThermoelectricMaterial object.

!syntax parameters /Materials/ThermoelectricMaterial

!syntax inputs /Materials/ThermoelectricMaterial

!syntax children /Materials/ThermoelectricMaterial
