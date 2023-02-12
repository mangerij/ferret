# Birefringence

!syntax description /AuxKernels/Birefringence

## Overview

Computes the difference between refractive indices $n_1(\mathbf{r})$ and $n_2(\mathbf{r})$ where $1$ and $2$ denote user-defined directions. Since this is an `AuxKernel`, these quantites are computed locally and can modulate in the simulation box if there are inhomogeneous fields (i.e. electric or elastic) that cause variations.

## Example Input File Syntax

!! Describe and include an example of how to use the Birefringence object.

!syntax parameters /AuxKernels/Birefringence

!syntax inputs /AuxKernels/Birefringence

!syntax children /AuxKernels/Birefringence
