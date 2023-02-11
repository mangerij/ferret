# FluctuationsIC

!alert construction title=Undocumented Class
The FluctuationsIC has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /ICs/FluctuationsIC

## Overview

Useful `ICs` function for benchmarking purposes. Typically in numerical simulations, the randomness of a function is prescribed by a random number generator which can depend on the processors. This IC allows for a pseudo-random function which will be the same on any machine.

\begin{equation}
  \begin{aligned}
 \delta_\eta = \cos{q_1 x - 4}\sin{q_1 y} + \cos{q_2 x} \cos{q_2 y} + ...
  \end{aligned}
\end{equation}



## Example Input File Syntax

!! Describe and include an example of how to use the FluctuationsIC object.

!syntax parameters /ICs/FluctuationsIC

!syntax inputs /ICs/FluctuationsIC

!syntax children /ICs/FluctuationsIC
