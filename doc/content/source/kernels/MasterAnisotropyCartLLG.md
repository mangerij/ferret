# MasterAnisotropyCartLLG

!alert construction title=Undocumented Class
The MasterAnisotropyCartLLG has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Kernels/MasterAnisotropyCartLLG

## Overview

Computes the residual and jacobian entries from the LLG equation due to a free energy density, $f_\mathrm{aniso}$, corresponding to magnetocrystalline anisotropy,

\begin{equation}
  \begin{aligned}
    f_\mathrm{aniso} = K_1 M_s^2 \left(\mathbf{m} \cdot \mathbf{w}\right)^2
  \end{aligned}
\end{equation}

where $\mathbf{m}$ is the normalized magnetization, $M_s$ is the saturation magnetization density, $K_1$ is the anisotropy constant ($>0$ for uniaxial, $<0$ for easy-plane), and $\mathbf{w}$ is the anisotropy director. Note that in the literature, it is typical to use $\mathbf{n}$ but we reserve this within the MOOSE ecosystem for the surface normals of the finite element mesh.

The effective field due to this term can be calculated, in index notation,

\begin{equation}
  \begin{aligned}
    H_\mathrm{aniso,k} &= - \frac{1}{\mu_0 M_s} \frac{\delta f_\mathrm{aniso}}{\delta m_k}\\
    &= - \frac{2 M_s}{\mu_0} w_k (m_j w_j)^2 \\
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

This equation for the $k^\mathrm{th}$ component in index notation is,

\begin{equation}
  \begin{aligned}
    \gamma \left(\psi_h, \epsilon_{ijk} m_i H_j \right) + \frac{\gamma \alpha}{1+\alpha^2}\left(\psi_h, \epsilon_{ijk} \epsilon_{rsj} m_i m_r H_s \right) = 0
  \end{aligned}
\end{equation}

with $\epsilon_{ijk}$ the Levi-Civita symbol. Inserting the expression for the effective field due to the magnetocrystalline anisotropy, we have the residual contribution for $m_k$,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_k} = 2 \gamma_\mathrm{H} K_1 \left(\psi_h, \epsilon_{ijk} m_i w_j (m_l w_l)^2 \right) + 2 \frac{\gamma_\mathrm{H} \alpha K_1}{1+\alpha^2} \left(\psi_h, \epsilon_{ijk} \epsilon_{rsj} m_i m_r w_s (m_l w_l)^2 \right).
  \end{aligned}
\end{equation}







## Example Input File Syntax

!! Describe and include an example of how to use the MasterAnisotropyCartLLG object.

!syntax parameters /Kernels/MasterAnisotropyCartLLG

!syntax inputs /Kernels/MasterAnisotropyCartLLG

!syntax children /Kernels/MasterAnisotropyCartLLG
