# RotoPolarCoupledEnergyPolarDerivativeAlt

!syntax description /Kernels/RotoPolarCoupledEnergyPolarDerivativeAlt

## Overview

The coupling free energy density between the electric polarization $\mathbf{P}$ and the oxygen octahedral antiphase tilt field $\mathbf{A}$ can be written, up to eighth order, as

\begin{equation}
  \begin{aligned}
    f_\mathrm{P-A} &= \tau_{ijkl} P_i P_j A_k A_l + \tau_{ijklmn}^a P_i P_j P_k A_m A_n + \tau_{ijklmn}^b P_i P_j A_k A_l A_m A_n \\
    &+\tau_{ijklmnpq}^c P_i P_j P_k P_l A_m A_n A_p A_q + \tau_{ijklmnpq}^d P_i P_j P_k P_l P_m P_n A_p A_q +\tau_{ijklmnpq}^e P_i P_j A_k A_l A_m A_n A_p A_q.
  \end{aligned}
\end{equation}

We need to compute variational derivatives of this expression with respect to the $\mathbf{P}$.

\begin{equation}
  \begin{aligned}
    f_\mathrm{P-A} &= 
  \end{aligned}
\end{equation}


## Example Input File Syntax

!! Describe and include an example of how to use the RotoPolarCoupledEnergyPolarDerivativeAlt object.

!syntax parameters /Kernels/RotoPolarCoupledEnergyPolarDerivativeAlt

!syntax inputs /Kernels/RotoPolarCoupledEnergyPolarDerivativeAlt

!syntax children /Kernels/RotoPolarCoupledEnergyPolarDerivativeAlt
