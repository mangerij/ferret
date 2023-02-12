# ScalarElasticEnergy

!syntax description /ScalarKernels/ScalarElasticEnergy

## Overview

Sets up the RHS of the coupled Landau-Khalatnikov equation. Since we are working with `ScalarKernels`, this is effectively a single grid point, or homogeneous calculation. In the context of micromagnetics, this is the macrospin analog for ferroelectrics. This object concerns the elastic energy contribution from the spontaneous strains (in the case of $\mathrm{BiFeO}_3$),

\begin{equation}
  \begin{aligned}
     \frac{\partial \varepsilon_{ij}}{\partial t} = - \Gamma_u \frac{\delta f_\mathrm{elastic}}{\delta \varepsilon_{ij}}
  \end{aligned}
\end{equation}

Ignoring the LHS and working in index notation,

\begin{equation}
  \begin{aligned}
     - \Gamma_u \frac{\delta f_\mathrm{elastic}}{\delta \varepsilon_{ij}} &= - \Gamma_u \frac{1}{2} \frac{\partial}{\partial \varepsilon_{ij}} \left( C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}\right) \\
  \end{aligned}
\end{equation}

In Voight notation, many of these terms vanish due to symmetry of the parent phase. The tensor $C_{ijkl}$ is the elastic stiffness matrix. Note that as opposed to a spatially dependent calculation, the strain tensor is assumed to be constant in space. This leads to our approximation that $\varepsilon_{ij}$ is a variable in this calculation as opposed to a postprocessed quantity depending on the elastic displacement vector $\mathbf{u}$.

## Example Input File Syntax

!! Describe and include an example of how to use the ScalarElasticEnergy object.

!syntax parameters /ScalarKernels/ScalarElasticEnergy

!syntax inputs /ScalarKernels/ScalarElasticEnergy

!syntax children /ScalarKernels/ScalarElasticEnergy
