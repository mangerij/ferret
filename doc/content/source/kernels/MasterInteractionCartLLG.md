# MasterInteractionCartLLG

!syntax description /Kernels/MasterInteractionCartLLG

## Overview

We consider the residual and jacobian contributions due to the demagnetization field (can also be an externally applied field). The free energy density due to the Zeeman interaction is,

\begin{equation}
  \begin{aligned}
    f_\mathrm{Zeeman} &= \mathbf{M} \cdot \mathbf{H} \\
    &= - M_s \mathbf{m} \cdot \nabla \Phi_\mathrm{H}.
  \end{aligned}
\end{equation}

This yields an effective field of

\begin{equation}
  \begin{aligned}
    \mathbf{H}_\mathrm{Zeeman} &= - \frac{1}{\mu_0 M_s} \frac{\delta f_\mathrm{Zeeman}}{\delta \mathbf{m}} \\
    &= \frac{1}{\mu_0} \nabla \Phi_\mathrm{H}. \\
  \end{aligned}
\end{equation}

The LLG equation for the normalized magnetization is,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{m}}{\partial t} = -\gamma \mathbf{m}\times \mathbf{H} - \frac{\gamma \alpha}{1+\alpha^2} \mathbf{m} \times \mathbf{m}\times \mathbf{H},
  \end{aligned}
\end{equation}

with $\gamma$ the gyromagnetic ratio. Multiplying by a test function $\psi_h$, moving over the RHS, neglecting the time derivative, and integrating over the volume, we have,

\begin{equation}
  \begin{aligned}
    \gamma\left(\psi_h, \mathbf{m}\times \mathbf{H}\right) + \frac{\gamma \alpha}{1+\alpha^2} \left(\psi_h, \mathbf{m} \times \mathbf{m}\times \mathbf{H}\right) = 0
  \end{aligned}
\end{equation}

This equation in index notation is,

\begin{equation}
  \begin{aligned}
    \gamma \left(\psi_h, \epsilon_{ijk} m_i H_j \right) + \frac{\gamma \alpha}{1+\alpha^2}\left(\psi_h, \epsilon_{ijk} \epsilon_{rsj} m_i m_r H_s \right) = 0
  \end{aligned}
\end{equation}

with $\epsilon_{ijk}$ the Levi-Civita symbol. Inserting the expression for the effective field due to the exchange stiffness, we have

\begin{equation}
  \begin{aligned}
    \gamma_\mathrm{H} \left(\psi_h, \epsilon_{ijk} m_i \frac{\partial \Phi_\mathrm{H}}{\partial x_j} \right) + \frac{\gamma_\mathrm{H} \alpha}{1+\alpha^2}\left(\psi_h, \epsilon_{ijk} \epsilon_{rsj} m_i m_r \frac{\partial \Phi_\mathrm{H}}{\partial x_s} \right) = 0
  \end{aligned}
\end{equation}

Therefore the residual for $m_k$ can be written as,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_k} = \gamma_\mathrm{H} \left(\psi_h, \epsilon_{ijk} m_i \frac{\partial \Phi_\mathrm{H}}{\partial x_j} \right) + \frac{\gamma_\mathrm{H} \alpha}{1+\alpha^2}\left(\psi_h, \epsilon_{ijk} \epsilon_{rsj} m_i m_r \frac{\partial \Phi_\mathrm{H}}{\partial x_s} \right).
  \end{aligned}
\end{equation}

The on-diagonal and off-diagonal jacobian expression is,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{m_k, m_\beta} &= \frac{\partial\mathcal{R}_{m_k}}{\partial m_\beta} \\
    &= \gamma_\mathrm{H} \left(\psi_h, \phi \epsilon_{\beta jk} \frac{\partial \Phi_\mathrm{H}}{\partial x_j} \right) + \frac{\gamma_\mathrm{H} \alpha}{1+\alpha^2}\left(\psi_h, \phi \epsilon_{\beta jk} \epsilon_{rsj} m_r \frac{\partial \Phi_\mathrm{H}}{\partial x_s} \right) + \frac{\gamma_\mathrm{H} \alpha}{1+\alpha^2}\left(\psi_h, \phi \epsilon_{ijk} \epsilon_{\beta sj} m_i \frac{\partial \Phi_\mathrm{H}}{\partial x_s} \right).
  \end{aligned}
\end{equation}

and

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{m_k, \Phi_\mathrm{H}} &= \frac{\partial\mathcal{R}_{m_k}}{\partial \Phi_\mathrm{H}} \\
    &= \gamma_\mathrm{H} \left(\psi_h, \epsilon_{ijk} m_i \frac{\partial \phi}{\partial x_j} \right) + \frac{\gamma_\mathrm{H} \alpha}{1+\alpha^2}\left(\psi_h, \epsilon_{ijk} \epsilon_{rsj} m_i m_r \frac{\partial \phi}{\partial x_s} \right) + \frac{\gamma_\mathrm{H} \alpha}{1+\alpha^2}\left(\psi_h, \epsilon_{ijk} \epsilon_{\beta sj} m_i m_r \frac{\partial \phi}{\partial x_s} \right).
  \end{aligned}
\end{equation}

where $\phi$ is the shape function of the finite element method.


## Example Input File Syntax

!! Describe and include an example of how to use the MasterInteractionCartLLG object.

!syntax parameters /Kernels/MasterInteractionCartLLG

!syntax inputs /Kernels/MasterInteractionCartLLG

!syntax children /Kernels/MasterInteractionCartLLG
