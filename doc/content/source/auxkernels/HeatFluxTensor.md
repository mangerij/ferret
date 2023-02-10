# HeatFluxTensor

!alert construction title=Undocumented Class
The HeatFluxTensor has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/HeatFluxTensor

## Overview

Computes the $m^\mathrm{th}$ component of the heat flux vector $\mathbf{q}$ (dependent on tensor quantities),

\begin{equation}
  \begin{aligned}
    q_m = \left(\Sigma_{i,j} \kappa_{ij} \frac{\partial T}{\partial x_m} + \Sigma_{i,j,k} S_{ij} \gamma_{jk} T \frac{\partial \Phi_\mathrm{E}}{\partial x_m} + \Sigma_{i,j,k,l} S_{ij} S_{jk} \gamma_{kl} T \frac{\partial T}{\partial x_m} \right)
  \end{aligned}
\end{equation}

The tensor components $\kappa_{ij}, S_{ij}$, and $\gamma_{ij}$ correspond to the thermal conductivity, Seebeck, and electrical conductivity respectively which require information of the material crystallographic directions.

## Example Input File Syntax

!! Describe and include an example of how to use the HeatFluxTensor object.

!syntax parameters /AuxKernels/HeatFluxTensor

!syntax inputs /AuxKernels/HeatFluxTensor

!syntax children /AuxKernels/HeatFluxTensor
