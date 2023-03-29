# PolarElectricPStrongEConst

!syntax description /Kernels/PolarElectricPStrongEConst

## Overview

Calculates the residual and jacobian terms related to the polarization vector $\mathbf{P}$ coupling to the electric field $\mathbf{E}$ in the time-dependent Landau-Ginzburg-Devosnhire (TDLGD) equation. The vector value of $\mathbf{E}$ is set in the input file and is not a MOOSE `Variable`. This is the distinction compared to the other `Kernel` named [PolarElectricPStrong](source/kernels/PolarElectricPStrong.md). We first consider the free energy density contribution,

\begin{equation}
  \begin{aligned}
    f_\mathrm{PE} &= - \mathbf{P} \cdot \mathbf{E} \\
  \end{aligned}
\end{equation}

where the TDLGD equation is,

\begin{equation}
  \begin{aligned}
    \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta f}{\delta \mathbf{P}},
  \end{aligned}
\end{equation}

The variational derivative of $f_\mathrm{PE}$ is, in index notation,

\begin{equation}
  \begin{aligned}
    \frac{\delta f_\mathrm{PE}}{\delta P_\eta} &= \frac{\partial f_\mathrm{PE}}{\partial P_\eta} \\
  \end{aligned}
\end{equation}

Ignoring the time derivative, bringing over the RHS to the LHS, and multiplying both sides by a test function $\psi_h$, we have the residual contribution for $P_\eta$,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{P_\eta} = \left(\psi_h, \frac{\partial f_\mathrm{PE}}{\partial P_\eta}\right) = \left(\psi_h, E_\eta\right).
  \end{aligned}
\end{equation}

The on- and off-diagonal jacobian entries,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{P_\eta,P_\beta} = \frac{\partial \mathcal{R}_{P_\eta}}{\partial P_\beta} = 0,
  \end{aligned}
\end{equation}

are zero since $\mathcal{R}_{P_\eta}$ does not explicitly depend on $\mathbf{P}$. In this alternative version of the [PolarElectricPStrong](source/kernels/PolarElectricPStrong.md) `Kernel`, $\mathbf{E}$ is just a vector value defined in the input file, so off-diagonal entries are also zero.

## Example Input File Syntax

!! Describe and include an example of how to use the PolarElectricPStrongEConst object.

!syntax parameters /Kernels/PolarElectricPStrongEConst

!syntax inputs /Kernels/PolarElectricPStrongEConst

!syntax children /Kernels/PolarElectricPStrongEConst
