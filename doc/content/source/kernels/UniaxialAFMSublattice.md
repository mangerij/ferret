# UniaxialAFMSublattice

!syntax description /Kernels/UniaxialAFMSublattice

## Overview

Here we consider the first example from the work of [!cite](Rezende2019) of a uniaxial antiferromagnet $\mathrm{MnF}_2$. The total free energy density of the material is,

\begin{equation}
  \begin{aligned}
    f = - H_0 \left(M_{1z} + M_{2z}\right) + \frac{H_\mathrm{E}}{M_s} \mathbf{M}_1 \cdot \mathbf{M}_2 - \frac{H_\mathrm{A}}{2 M_s} \left(M_{1z}^2 + M_{2z}^2\right).
  \end{aligned}
\end{equation}

We look for normalized solutions for $\mathbf{m}_\eta = \mathbf{M}_\eta / M_s$ to the two-sublattice LLG-LLB equation in the presence of zero Gilbert damping ($\alpha = 0$),

\begin{equation}\label{eqn:LLG_LLB}
  \begin{aligned}
    \frac{d \mathbf{m}_\eta}{dt} = -\gamma \left(\mathbf{m}_\eta \times \mathbf{H}_\eta\right) + \frac{\gamma\tilde\alpha_\parallel}{(1+\alpha^2)} m_{\eta}^2 \left[ m_{\eta}^2 - 1 \right]\mathbf{m}_\eta,
  \end{aligned}
\end{equation}

We ignore the time dependence and the longitudinal damping term (from the LLB approximation). Then, we multiply both sides by a test function $\psi_h$ and integrate over the volume, which yields,

\begin{equation}
  \begin{aligned}
    -\left(\psi_h, \gamma \mathbf{m}_\eta \times \mathbf{H}_\eta \right) = 0
  \end{aligned}
\end{equation}

The effective field is computed as,

\begin{equation}
  \begin{aligned}
    \mathbf{H}_\eta = -\frac{1}{\mu_0 M_s}\left(\frac{\delta f}{\delta \mathbf{m}_\eta}\right).
  \end{aligned}
\end{equation}

This gives the residual contribution (in index notation) for $m_{\eta k}$ to be

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_{\eta k}} &= -\gamma \left(\psi_h, \epsilon_{ijk} m_{i\eta} H_{j\eta} \right)\\
    &= \frac{\gamma}{\mu_0 M_s} \left(\psi_h, \epsilon_{ijk} m_{i\eta} \frac{\delta f}{\delta m_{\eta j}} \right)\\
  \end{aligned}
\end{equation}

The variational derivatives of $f$ can be computed as,

\begin{equation}
  \begin{aligned}
    \frac{\delta f}{\delta m_{\eta j}} &= M_s \frac{\partial}{\partial m_{\eta j}} \left(- H_0 \left(m_{1z} + m_{2z}\right) + H_\mathrm{E} m_{1 \beta} m_{2 \beta} - \frac{H_\mathrm{A}}{2} \left(m_{1z}^2 + m_{2z}^2\right)\right)\\
    &= M_s \left(- H_0 \delta_{zj} + H_\mathrm{E} \delta_{\beta j} m_{\eta \beta} - H_\mathrm{A} \delta_{zj} m_{\eta z}\right)\\
  \end{aligned}
\end{equation}

Therefore,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_{\eta k}} &= -\frac{\gamma}{\mu_0} \left(\psi_h, \epsilon_{ijk} m_{i\eta} \left\{- H_0 \delta_{zj} + H_\mathrm{E} \delta_{\beta j} m_{\eta \beta} - H_\mathrm{A} \delta_{zj} m_{\eta z}\right\} \right)\\
  \end{aligned}
\end{equation}

## Example Input File Syntax

!! Describe and include an example of how to use the UniaxialAFMSublattice object.

!syntax parameters /Kernels/UniaxialAFMSublattice

!syntax inputs /Kernels/UniaxialAFMSublattice

!syntax children /Kernels/UniaxialAFMSublattice
