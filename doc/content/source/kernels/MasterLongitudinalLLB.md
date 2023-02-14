# MasterLongitudinalLLB

!syntax description /Kernels/MasterLongitudinalLLB

## Overview

We consider the residual and jacobian contributions due to the restoring force from the LLG-LLB approximation. The LLG-LLB equation is,

\begin{equation}\label{LLG-LLB}
  \begin{aligned}
    \frac{\partial \mathbf{M}}{\partial t} = -\gamma \mathbf{M}\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M} \times \mathbf{M}\times \mathbf{H} - \frac{\gamma \alpha_\mathrm{LLB}}{1+\alpha^2} m^2 \left(m^2 - 1\right) \mathbf{m}
  \end{aligned}
\end{equation}

where we write in normalized form,

\begin{equation}
  \begin{aligned}
    \frac{\partial \mathbf{m}}{\partial t} = -\gamma \mathbf{m}\times \mathbf{H} - \frac{\gamma \alpha}{1+\alpha^2} \mathbf{m} \times \mathbf{m}\times \mathbf{H} - \frac{\gamma \tilde{\alpha}_\mathrm{LLB}}{1+\alpha^2} m^2 \left(m^2 - 1\right) \mathbf{m}
  \end{aligned}
\end{equation}

where $\tilde{\alpha}_\mathrm{LLB} = \alpha_\mathrm{LLB} / M_s$ is the longitudinal damping constant. Moving over the RHS, multiplying both sides by a test function $\psi_h$, integrating over the volume, and isolating the the LLB term we have,

\begin{equation}
  \begin{aligned}
    \frac{\gamma \tilde{\alpha}_\mathrm{LLB}}{1+\alpha^2} \left(\psi_h, m^2 \left(m^2 - 1\right) \mathbf{m} \right) = 0.
  \end{aligned}
\end{equation}

Therefore, in index notation, the residual contribution for $m_\eta$ is simply,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_\eta} = \frac{\gamma \tilde{\alpha}_\mathrm{LLB}}{1+\alpha^2} \left(\psi_h, m_j m_j \left(m_k m_k - 1\right) m_\eta \right).
  \end{aligned}
\end{equation}

with the on- and off-diagonal jacobian contributions calculated via,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{m_\eta, m_\beta} &= \frac{\partial \mathcal{R}_{m_\eta}}{\partial m_\beta} \\
    &= \frac{\gamma \tilde{\alpha}_\mathrm{LLB}}{1+\alpha^2} \frac{\partial}{\partial m_\beta}\left(\psi_h, m_j m_j\left(m_k m_k - 1\right) m_\eta \right) \\
    &= \frac{\gamma \tilde{\alpha}_\mathrm{LLB}}{1+\alpha^2} \left(\psi_h, 4 \phi m_j m_j m_\eta m_\beta - 2 \phi m_\beta m_\eta - \phi m_k m_k \delta_{\eta \beta} \right) \\
  \end{aligned}
\end{equation}

with $\phi$ the shape function from the finite element method. Note that the last term of this equation is only nonzero if we are computing the on-diagonal jacobian entry in the case of $\beta = \eta$.

## Example Input File Syntax

!! Describe and include an example of how to use the MasterLongitudinalLLB object.

!syntax parameters /Kernels/MasterLongitudinalLLB

!syntax inputs /Kernels/MasterLongitudinalLLB

!syntax children /Kernels/MasterLongitudinalLLB
