# HarmonicFieldAux

!alert construction title=Undocumented Class
The HarmonicFieldAux has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/HarmonicFieldAux

## Overview

Computes a harmonic field of the form,

\begin{equation}
  \begin{aligned}
    A = A_0 c_f \sin{\left(\omega \left(t-t_\mathrm{shift}\right)}
  \end{aligned}
\end{equation}

where $A_0$ and $c_f$ are amplitudes and correction factors (i.e. 1/\sqrt{2}). The frequency is $\omega$, simulation time $t$, and temporal phase shift (if needed) is $t_\mathrm{shift}$. Note that it is also possible to utilized this object to turn the field on at a time $t_\mathrm{on}$ and off at a time $t_\mathrm{off}$.

## Example Input File Syntax

!! Describe and include an example of how to use the HarmonicFieldAux object.

!syntax parameters /AuxKernels/HarmonicFieldAux

!syntax inputs /AuxKernels/HarmonicFieldAux

!syntax children /AuxKernels/HarmonicFieldAux
