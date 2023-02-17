# ElasticEnergyDensity

!syntax description /AuxKernels/ElasticEnergyDensity

## Overview

Computes the elastic energy density due to

\begin{equation}
  \begin{aligned}
    f_\mathrm{elastic} &= \frac{1}{2} \sigma_{kl} \varepsilon_{kl} \\
    &= \frac{1}{2} C_{ijkl} \varepsilon_{ij} \varepsilon_{kl}.
  \end{aligned}
\end{equation}

where $\sigma_{kl}$ and $\varepsilon_{kl}$ is the elastic stress and strain tensor component respectively. Shown here is an equivalent expression involving the elastic stiffness tensor $C_{ijkl}$ of rank four.

## Example Input File Syntax

!! Describe and include an example of how to use the ElasticEnergyDensity object.

!syntax parameters /AuxKernels/ElasticEnergyDensity

!syntax inputs /AuxKernels/ElasticEnergyDensity

!syntax children /AuxKernels/ElasticEnergyDensity
