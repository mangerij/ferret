# HarmonicFieldAux

!syntax description /AuxKernels/HarmonicFieldAux

## Overview

Computes a harmonic field of the form,

\begin{equation}
  \begin{aligned}
    A = A_0 c_f \sin{\omega (t - t_\mathrm{shift})}
  \end{aligned}
\end{equation}

where $A_0$ and $c_f$ are amplitudes and correction factors (i.e. $1/\sqrt{2}$). The frequency is $\omega$, simulation time $t$, and temporal phase shift (if needed) is $t_\mathrm{shift}$. Note that it is also possible to utilized this object to turn the field on at a time $t_\mathrm{on}$ and off at a time $t_\mathrm{off}$.

## Example Input File Syntax

!! Describe and include an example of how to use the HarmonicFieldAux object.

!syntax parameters /AuxKernels/HarmonicFieldAux

!syntax inputs /AuxKernels/HarmonicFieldAux

!syntax children /AuxKernels/HarmonicFieldAux
