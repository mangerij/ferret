# GlobalBFOMaterialRVEUserObject

!syntax description /UserObjects/GlobalBFOMaterialRVEUserObject

## Overview

Computes the necessary eigenstress needed for the `GlobalStrain` system which performs the calculation,

\begin{equation}
  \begin{aligned}
    \int\limits_\Omega \sigma_{ij}^\mathrm{total} d^3 \mathbf{r} = 0.
  \end{aligned}
\end{equation}

over the computational volume $\Omega$. For ferroelectric $\mathrm{BiFeO}_3$, the eigenstrain that is generated via the electrostrictive *and* rotostrictive coupling is,

\begin{equation}
  \begin{aligned}
    \varepsilon_{ij}^\mathrm{eig} = Q_{ijkl} P_k P_l + R_{ijkl} A_k A_l.
  \end{aligned}
\end{equation}

where $P_k$ and $P_l$ are components of the polarization vector $\mathbf{P}$ and $A_k$ and $A_l$ are components of the oxygen octahedral antiphase tilt vector $\mathbf{A}$. The quantities $Q_{ijkl}$ and $R_{ijkl}$ are the electrostrictive and rotostrictive tensors of rank four. Therefore the eigenstress is,

\begin{equation}
  \begin{aligned}
    \sigma_{ij}^\mathrm{eig} = C_{ijkl} \left(Q_{klmn} P_m P_n + R_{klmn} A_m A_n\right).
  \end{aligned}
\end{equation}

The total stress $\sigma_{ij}^\mathrm{total}$ is,

\begin{equation}
  \begin{aligned}
    \sigma_{ij}^\mathrm{total} = \sigma_{ij} + \sigma_{ij}^\mathrm{eig}.
  \end{aligned}
\end{equation}

More details of this approach are provided in [!cite](Biswas2020).

## Example Input File Syntax

!! Describe and include an example of how to use the GlobalBFOMaterialRVEUserObject object.

!syntax parameters /UserObjects/GlobalBFOMaterialRVEUserObject

!syntax inputs /UserObjects/GlobalBFOMaterialRVEUserObject

!syntax children /UserObjects/GlobalBFOMaterialRVEUserObject
