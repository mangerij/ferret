# PolarElectricEStrong

!syntax description /Kernels/PolarElectricEStrong

## Overview

Computes the residual and jacobian terms due to the bound polarization charge $\rho_b = \nabla \cdot \mathbf{P}$ (RHS of the Poisson equation). Consider the Poisson equation,

\begin{equation}
  \begin{aligned}
    - \frac{\partial}{\partial x_j} \epsilon_b \frac{\partial \Phi_\mathrm{E}}{\partial x_j} = \frac{\partial P_k}{\partial x_k}
  \end{aligned}
\end{equation}

where $\mathbf{P}$ is the electric polarization vector, $\Phi_\mathrm{E}$ is the electrostatic potential, and $\epsilon_b$ is a background relative permittivity. Multiplying both sides by a test function suitable for finite element analysis, $\psi_h$ and integrating over both sides, we have,

\begin{equation}
  \begin{aligned}
    - \left(\psi_h, \frac{\partial}{\partial x_j} \epsilon_b \frac{\partial \Phi_\mathrm{E}}{\partial x_j}\right) - \left(\psi_h, \frac{\partial P_k}{\partial x_k}\right) = 0
  \end{aligned}
\end{equation}

Integrating by parts,

\begin{equation}
  \begin{aligned}
    \left(\frac{\partial \psi_h}{\partial x_j}, \epsilon_b \frac{\partial \Phi_\mathrm{E}}{\partial x_j}\right) + \left\langle\frac{\partial \psi_h}{\partial x_j}, \epsilon_b \frac{\partial \Phi_\mathrm{E}}{\partial x_j}\right\rangle - \left(\frac{\psi_h}{\partial x_k}, P_k\right) + \left\langle\frac{\psi_h}{\partial x_k}, P_k\right\rangle= 0
  \end{aligned}
\end{equation}

with surface terms $\langle , \rangle$ vanish due to properties of $\psi_h$. Ignoring the Laplace term, which is handled by the `Electrostatics` `Kernel`, we have the residual contribution for $\Phi_\mathrm{E}$,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{\Phi_\mathrm{E}} = -\left(\frac{\psi_h}{\partial x_k}, P_k\right).
  \end{aligned}
\end{equation}

This does not explicitly depend on $\Phi_\mathrm{E}$ so therefore on-diagonal jacobians,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{\Phi_\mathrm{E},\Phi_\mathrm{E}} = \frac{\partial \mathcal{R}_{\Phi_\mathrm{E}}}{\partial \Phi_\mathrm{E}} = 0,
  \end{aligned}
\end{equation}

are zero. We need to compute the off-diagonal jacobian entries due to,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{\Phi_\mathrm{E}, P_\eta} &= -\frac{\partial \mathcal{R}_{\Phi_\mathrm{E}}}{\partial P_\eta}\\
    &= \left(\frac{\psi_h}{\partial x_k}, \delta_{k\eta} \phi \right),
  \end{aligned}
\end{equation}

since $\partial P_k / \partial P_\eta = \delta_{k\eta} \phi$ with $\phi$ a shape function of the finite element method.

## Example Input File Syntax

!! Describe and include an example of how to use the PolarElectricEStrong object.

!syntax parameters /Kernels/PolarElectricEStrong

!syntax inputs /Kernels/PolarElectricEStrong

!syntax children /Kernels/PolarElectricEStrong
