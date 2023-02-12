# DomainVariantPopulation

!syntax description /Postprocessors/DomainVariantPopulation

## Overview

Computes the fraction of the domain variation population in the compuational box of volume $\Omega$,

\begin{equation}
  \begin{aligned}
     \xi_j = \frac{1}{\Omega}\int\limits_\Omega d^3\mathbf{r} \left\{ \frac{|P_j|}{\sqrt{P_x^2 + P_y^2 + P_z^2}} \right\}.
  \end{aligned}
\end{equation}

This expression is only valid for tetragonal ferroelectrics but we plan to generalize this to orthorhombic and rhombohedral symmetries in the future.

## Example Input File Syntax

!! Describe and include an example of how to use the DomainVariantPopulation object.

!syntax parameters /Postprocessors/DomainVariantPopulation

!syntax inputs /Postprocessors/DomainVariantPopulation

!syntax children /Postprocessors/DomainVariantPopulation
