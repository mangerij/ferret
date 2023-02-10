# ElasticEnergyDensity

!alert construction title=Undocumented Class
The ElasticEnergyDensity has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/ElasticEnergyDensity

## Overview

Computes the elastic energy density locally due to

\begin{equation}
  \begin{aligned}
    f_\mathrm{elastic} = \frac{1}{2} \sigma_{kl} \varepsilon_{kl} = \frac{1}{2} C_{ijkl} \varepsilon_{ij} \varepsilon_{kl}.
  \end{aligned}
\end{equation}

## Example Input File Syntax

!! Describe and include an example of how to use the ElasticEnergyDensity object.

!syntax parameters /AuxKernels/ElasticEnergyDensity

!syntax inputs /AuxKernels/ElasticEnergyDensity

!syntax children /AuxKernels/ElasticEnergyDensity
