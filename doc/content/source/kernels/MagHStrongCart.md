# MagHStrongCart

!syntax description /Kernels/MagHStrongCart

## Overview

The magnetostatic Poisson equation reads,

\begin{equation}
  \begin{aligned}
    \frac{\partial^2\Phi_\mathrm{H}}{\partial x_j^2} &= \frac{\partial M_k}{\partial x_k}, \\
    &= M_s\frac{\partial m_k}{\partial x_k}
  \end{aligned}
\end{equation}

where $\Phi_\mathrm{H}$ is the magnetostatic potential due to a bound magnetization charge $\nabla \cdot \mathbf{M}$. Multiplying both sides by a test function $\psi_h$ and integrating over the computational volume, we have

\begin{equation}
  \begin{aligned}
    \left(\psi_h,\frac{\partial^2\Phi_\mathrm{H}}{\partial x_j^2}\right) - \left(M_s \psi_h \frac{\partial m_k}{\partial x_k}\right) = 0.
  \end{aligned}
\end{equation}

Integrating by parts,

\begin{equation}
  \begin{aligned}
    -\left(\frac{\partial\psi_h}{\partial x_j},\frac{\partial\Phi_\mathrm{H}}{\partial x_j}\right) + \left\langle\frac{\partial\psi_h}{\partial x_j},\frac{\partial\Phi_\mathrm{H}}{\partial x_j}\right\rangle  - \left(M_s \frac{\partial\psi_h}{\partial x_k}, m_k\right) + \left\langle M_s \frac{\partial\psi_h}{\partial x_k}, m_k\right\rangle= 0.
  \end{aligned}
\end{equation}

The surface terms $\langle m \rangle$ vanish due to properties of $\psi_h$ and we are left with the residual contribution for $\Phi_\mathrm{H}$.

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{\Phi_\mathrm{H}} = \underbrace{-\left(\frac{\partial\psi_h}{\partial x_j},\frac{\partial\Phi_\mathrm{H}}{\partial x_j}\right)}_\mathrm{Electrostatics} \,\,\, \underbrace{- \left(M_s \frac{\partial\psi_h}{\partial x_k}, m_k\right)}_\mathrm{MagHStrongCart}.
  \end{aligned}
\end{equation}

Note that the first term is just the Laplacian which we name `Electrostatics` `Kernel` with $\epsilon_b = 1$. The on-diagonal jacobian entries computed by,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{\Phi_\mathrm{H},\Phi_\mathrm{H}} = \frac{\partial \mathcal{R}_{\Phi_\mathrm{H}}}{\partial \Phi_\mathrm{H}} = 0
  \end{aligned}
\end{equation}

are zero since this term does not explicitly depend on $\Phi_\mathrm{H}$. However, we still need to compute the off-diagonal jacobian entries due to variations with respect to $m_\eta$

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{\Phi_\mathrm{H},m_\eta} &= \frac{\partial \mathcal{R}_{\Phi_\mathrm{H}}}{\partial m_\eta} \\
    &= \frac{\partial}{\partial m_\eta} \left(M_s \frac{\partial\psi_h}{\partial x_k}, m_k\right)\\
    &= \left(M_s \frac{\partial\psi_h}{\partial x_k}, \delta_{k\eta}\phi\right)
  \end{aligned}
\end{equation}

where $\phi$ is a shape function of the finite element method and $\delta_{k\eta}$ the Kronecker product.

## Example Input File Syntax

!! Describe and include an example of how to use the MagHStrongCart object.

!syntax parameters /Kernels/MagHStrongCart

!syntax inputs /Kernels/MagHStrongCart

!syntax children /Kernels/MagHStrongCart
