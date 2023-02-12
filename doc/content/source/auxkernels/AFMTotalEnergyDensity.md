# AFMTotalEnergyDensity

!syntax description /AuxKernels/AFMTotalEnergyDensity

## Overview

Computes the sum of the energy densities relevant for the antiferromagnetic (AFM) ferroelectric (i.e. $\mathrm{BiFeO}_3$),

\begin{equation}
  \begin{aligned}
    f_\mathrm{total} &=  f_\mathrm{super} + f_\mathrm{DMI} + f_\mathrm{easy} + f_\mathrm{anis} + f_\mathrm{exch}
  \end{aligned}
\end{equation}

where $f_\mathrm{super}$, $f_\mathrm{DMI}$, $f_\mathrm{easy}$, $f_\mathrm{anis}$, and $f_\mathrm{exch}$ are the short-range superexchange, Dzyaloshinskii–Moriya interaction, easy-plane anisotropy, weak in-plane anisotropy, and long-range AFM exchange.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMTotalEnergyDensity object.

!syntax parameters /AuxKernels/AFMTotalEnergyDensity

!syntax inputs /AuxKernels/AFMTotalEnergyDensity

!syntax children /AuxKernels/AFMTotalEnergyDensity
