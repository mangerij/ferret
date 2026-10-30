# MasterAnisotropyCartLLG

!syntax description /Kernels/MasterAnisotropyCartLLG

## Overview

Computes the residual and jacobian entries from the LLG equation due to a free energy density, $f_\mathrm{aniso}$, corresponding to uniaxial magnetocrystalline anisotropy,

\begin{equation}
  \begin{aligned}
    f_\mathrm{aniso} = - K_1 \left(\mathbf{m} \cdot \mathbf{w}\right)^2 - K_2 \left(\mathbf{m} \cdot \mathbf{w}\right)^4
  \end{aligned}
\end{equation}

where $\mathbf{m}$ is the normalized magnetization, $K_1$ and $K_2$ are the anisotropy constants and $\mathbf{w}$ is the (normalized) anisotropy director. This is the standard convention (as in OOMMF, MuMax3 and the $\mu$MAG standard problems): $K_1 > 0$ gives an easy axis along $\mathbf{w}$, $K_1 < 0$ an easy plane normal to $\mathbf{w}$, and $K_1 > 0$ with $K_2 < -K_1/2$ an easy cone at $\cos^2\theta = K_1/(2|K_2|)$. The material property `K2` is optional and is taken as zero when no material provides it. Note that in the literature, it is typical to use $\mathbf{n}$ but we reserve this within the MOOSE ecosystem for the surface normals of the finite element mesh.

!alert warning title=Sign convention changed
Before October 2026 this kernel implemented $f_\mathrm{aniso} = + K_1 (\mathbf{m}\cdot\mathbf{w})^2$, so input files written for earlier versions used $K_1 = -K_u$ for an easy axis. Such inputs must flip the sign of `K1`. The energy postprocessor [MasterMagneticAnisotropyEnergy](MasterMagneticAnisotropyEnergy.md) already used $-K_1$ and is consistent with this kernel, including the $K_2$ term.

The effective field due to this term is, in index notation,

\begin{equation}
  \begin{aligned}
    H_\mathrm{aniso,k} &= - \frac{1}{\mu_0 M_s} \frac{\partial f_\mathrm{aniso}}{\partial m_k}\\
    &= \frac{1}{\mu_0 M_s} \left[2 K_1 \left(m_j w_j\right) + 4 K_2 \left(m_j w_j\right)^3\right] w_k.
  \end{aligned}
\end{equation}

The LLG equation for the normalized magnetization, in the Landau-Lifshitz form used by all `Master*CartLLG` kernels, is

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{m}}{\partial t} = -\frac{\gamma}{1+\alpha^2} \left[\mathbf{m}\times \mathbf{H} + \alpha\, \mathbf{m} \times \left(\mathbf{m}\times \mathbf{H}\right)\right],
  \end{aligned}
\end{equation}

with $\gamma$ the gyromagnetic ratio (`g0`) and $\alpha$ the Gilbert damping. Multiplying by a test function $\psi_h$ and integrating over the volume, the kernel supplies the term

\begin{equation}
  \begin{aligned}
    \frac{\gamma}{1+\alpha^2}\left[\left(\psi_h, \epsilon_{kij} m_i H_j \right) + \alpha \left(\psi_h, \left(m_k m_j H_j - m_j m_j H_k\right) \right)\right],
  \end{aligned}
\end{equation}

where $\epsilon_{ijk}$ is the Levi-Civita symbol and $|\mathbf{m}|^2$ is not replaced by 1. Inserting the anisotropy field gives the residual contribution for $m_k$,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_k} = \frac{\gamma}{\mu_0 H_\mathrm{scale} M_s (1+\alpha^2)} \left(\psi_h, \left[2 K_1 \left(m_l w_l\right) + 4 K_2 \left(m_l w_l\right)^3\right] \left[\epsilon_{kij} m_i w_j + \alpha \left(m_k m_j w_j - m_j m_j w_k\right) \right] \right),
  \end{aligned}
\end{equation}

with $H_\mathrm{scale}$ the optional field scale (`Hscale`).

## Hard axis and the Landau-Lifshitz-Bloch route

The optional hard-axis term $K_h(\mathbf{m}\cdot\mathbf{n}_h)^2$ (material properties `Kh`, `nhx`, `nhy`, `nhz`; $K_h > 0$ makes $\mathbf{n}_h$ hard) works in both routes and gives biaxial walls with a Walker breakdown.

With `llb = true` the kernel solves the LLB equation ([Finite-temperature micromagnetics](theory/intro_LLB.md)) with the transverse form of the energy, $-K_1[(\mathbf{m}\cdot\mathbf{n})^2 - |\mathbf{m}|^2] - K_2[(\mathbf{m}\cdot\mathbf{n})^4 - |\mathbf{m}|^4]$. It has the same torque as the usual form but does not shift $|\mathbf{m}|$ away from $m_e$ in the domains. `K1` is $K(T)/m_e^2$; the material properties `alpha_perp` and `alpha_par` are required.

!listing test/tests/magnets/llb_biaxial_chi.i block=Kernels/MasterAnisotropyCartLLG_x

## Example Input File Syntax

!! Describe and include an example of how to use the MasterAnisotropyCartLLG object.

!syntax parameters /Kernels/MasterAnisotropyCartLLG

!syntax inputs /Kernels/MasterAnisotropyCartLLG

!syntax children /Kernels/MasterAnisotropyCartLLG
