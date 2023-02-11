# GlobalATiO3MaterialRVEUserObject

!alert construction title=Undocumented Class
The GlobalATiO3MaterialRVEUserObject has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /UserObjects/GlobalATiO3MaterialRVEUserObject

## Overview

Computes the necessary eigenstress needed for the `GlobalStrain` system which performs the calculation,

\begin{equation}
  \begin{aligned}
    \int\limits_\Omega \sigma_{ij}^\mathrm{total} d^3 \mathbf{r} = 0.
  \end{aligned}
\end{equation}

over the computational volume $\Omega$. For a conventional ferroelectric material (i.e. $\mathrm{ATiO}_3$ where A is a cation such as Pb or Ba), the eigenstrain that is generated via the electrostrictive coupling is,

\begin{equation}
  \begin{aligned}
    \varepsilon_{ij}^\mathrm{eig} = Q_{ijkl} P_k P_l.
  \end{aligned}
\end{equation}

where $P_k$ and $P_l$ are components of the polarization vector $\mathbf{P}$. The quantites $Q_{ijkl}$ is the components of the electrostrictive tensor of rank four. Therefore the eigenstress is,

\begin{equation}
  \begin{aligned}
    \sigma_{ij}^\mathrm{eig} = C_{ijkl} Q_{klmn} P_m P_n.
  \end{aligned}
\end{equation}

The total stress $\sigma_{ij}^\mathrm{total}$ is calculated via,

\begin{equation}
  \begin{aligned}
    \sigma_{ij}^\mathrm{total} = \sigma_{ij} + \sigma_{ij}^\mathrm{eig}.
  \end{aligned}
\end{equation}

## Example Input File Syntax

!! Describe and include an example of how to use the GlobalATiO3MaterialRVEUserObject object.

!syntax parameters /UserObjects/GlobalATiO3MaterialRVEUserObject

!syntax inputs /UserObjects/GlobalATiO3MaterialRVEUserObject

!syntax children /UserObjects/GlobalATiO3MaterialRVEUserObject
