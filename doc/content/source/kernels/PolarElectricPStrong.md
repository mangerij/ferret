# PolarElectricPStrong

!syntax description /Kernels/PolarElectricPStrong

## Overview

Calculates the residual and jacobian terms related to the polarization vector $\mathbf{P}$ coupling to the electric field $\mathbf{E} = - \nabla \Phi_\mathrm{E}$ in the time-dependent Landau-Ginzburg-Devosnhire (TDLGD) equation. We first consider the free energy density contribution,

\begin{equation}
  \begin{aligned}
    f_\mathrm{PE} &= - \mathbf{P} \cdot \mathbf{E} \\
    &= \mathbf{P} \cdot \nabla \Phi_\mathrm{E},
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
    \mathcal{R}_{P_\eta} = \left(\psi_h, \frac{\partial f_\mathrm{PE}}{\partial P_\eta}\right) = \left(\psi_h, \frac{\partial \Phi_\mathrm{E}}{\partial x_\eta}\right).
  \end{aligned}
\end{equation}

The on- and off-diagonal jacobian entries,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{P_\eta,P_\beta} = \frac{\partial \mathcal{R}_{P_\eta}}{\partial P_\beta} = 0,
  \end{aligned}
\end{equation}

are zero since $\mathcal{R}_{P_\eta}$ does not explicitly depend on $\mathbf{P}$. We, however, do need to calculate the off-diagonal entry corresponding to,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{P_\eta,\Phi_\mathrm{E}} &= \frac{\partial \mathcal{R}_{P_\eta}}{\partial \Phi_\mathrm{E}}\\
    &= \frac{\partial}{\partial \Phi_\mathrm{E}} \left(\psi_h, \frac{\partial \Phi_\mathrm{E}}{\partial x_\eta}\right)\\
    &= \left(\psi_h, \frac{\partial \phi}{\partial x_\eta}\right)\\
  \end{aligned}
\end{equation}

since $\partial^2 \Phi_\mathrm{E} / \partial x_\eta \partial \Phi_\mathrm{E} = \partial \phi / \partial x_\eta$ with $\phi$ the shape function of the finite element method.


## Example Input File Syntax

!! Describe and include an example of how to use the PolarElectricPStrong object.

!syntax parameters /Kernels/PolarElectricPStrong

!syntax inputs /Kernels/PolarElectricPStrong

!syntax children /Kernels/PolarElectricPStrong
