# ElasticEnergy

!syntax description /Postprocessors/ElasticEnergy

## Overview

Computes the elastic energy of the simulation box (volume $\Omega$) due to

\begin{equation}
  \begin{aligned}
    F_\mathrm{elastic} &= \int\limits_\Omega \left\{ \frac{1}{2} \sigma_{kl} \varepsilon_{kl} \right\}\\
    &= \int\limits_\Omega \left\{\frac{1}{2} C_{ijkl} \varepsilon_{ij} \varepsilon_{kl} \right\}.
  \end{aligned}
\end{equation}

where $\sigma_{kl}$ and $\varepsilon_{kl}$ is the elastic stress and strain tensor component respectively. Shown here is an equivalent expression involving the elastic stiffness tensor $C_{ijkl}$ of rank four.

## Example Input File Syntax

!! Describe and include an example of how to use the ElasticEnergy object.

!syntax parameters /Postprocessors/ElasticEnergy

!syntax inputs /Postprocessors/ElasticEnergy

!syntax children /Postprocessors/ElasticEnergy
