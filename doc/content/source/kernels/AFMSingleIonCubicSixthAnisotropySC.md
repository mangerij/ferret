# AFMSingleIonCubicSixthAnisotropySC

!syntax description /Kernels/AFMSingleIonCubicSixthAnisotropySC

## Overview

Consider the following sixth-order free energy contribution,

\begin{equation}
  \begin{aligned}
    f_\mathrm{aniso} = \sum\limits_{\eta = 1,2} \left(K_1^c + a |\mathbf{A}|\right) m_{x\eta}^2 m_{y\eta}^2 m_{z\eta}^2
  \end{aligned}
\end{equation}

where $\mathbf{m}_\eta$ is the normalized sublattice magnetization. The vector $\mathbf{A}$ is the oxygen octahedral antiphase tilt vector (for the $\mathrm{BiFeO}_3$ material). The suffix `SC` in this `Kernel` stands for strong-coupling where we treat $\mathbf{A}$ as a MOOSE `Variable`. Therefore, we will need to compute additional off-diagonal jacobian entries. The effective field due to this term is,

\begin{equation}
  \begin{aligned}
    \mathbf{H}^\mathrm{aniso}_\eta &= - \frac{1}{\mu_0 M_s} \frac{\delta f_\mathrm{aniso}}{\delta \mathbf{m}_\eta}.
  \end{aligned}
\end{equation}

The effective field is evaluated at every time step of the two sublattice LLG-LLB equation,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{m}}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \mathbf{m}\times \mathbf{H} - \frac{\gamma \alpha}{1+\alpha^2} \mathbf{m} \times \mathbf{m}\times \mathbf{H},
  \end{aligned}
\end{equation}

here $\gamma$ is the electron gyromagnetic factor and $\alpha$ the phenomenological Gilbert damping parameter. Ignoring the time derivative, moving over the RHS to the LHS, and multiplying by a test function $\psi_h$, we have

\begin{equation}
  \begin{aligned}
    \frac{\gamma}{(1+\alpha^2)}\left(\psi_h, \mathbf{m}_\eta \times \mathbf{H}_\eta\right) + \frac{\gamma \alpha}{1+\alpha^2} \left(\psi_h, \mathbf{m}_\eta \times \mathbf{m}_\eta \times \mathbf{H}_\eta \right) = 0
  \end{aligned}
\end{equation}

Substituting the expression for $\mathbf{H}_\eta^\mathrm{aniso}$ and writing the above in index notation gives the residual contribution for the $m_{k\eta}$ component,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_{k\eta}} = \frac{\gamma}{\mu_0 M_s (1+\alpha^2)}\left\{\left(\psi_h, \epsilon_{ijk} m_{i\eta} \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} \right) + \alpha \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} m_{i\eta} m_{j\eta} \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} \right)\right\}.
  \end{aligned}
\end{equation}

The on- and off-diagonal jacobian entries are given by,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{m_{k\eta},m_{\beta\eta}} &= \frac{\partial \mathcal{R}_{m_{k\eta}}}{\partial m_{\beta\eta}} = \frac{\gamma}{\mu_0 M_s (1+\alpha^2)} \frac{\partial}{\partial m_{\beta\eta}} \left\{\left(\psi_h, \epsilon_{ijk} m_{i\eta} \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} \right) + \alpha \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} m_{i\eta} m_{j\eta} \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} \right)\right\} \\
    &=  \frac{\gamma}{\mu_0 M_s (1+\alpha^2)} \left\{\left(\psi_h, \epsilon_{ijk} \left[ \delta_{i\beta} \phi \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} + m_{i\eta} \frac{\partial^2 f_\mathrm{aniso}}{\partial m_{k\eta}\partial m_{\beta\eta}} \right] \right) + \alpha \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} \left[ \delta_{i\beta} \phi m_{j\eta} \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} + m_{i\eta} \delta_{j\beta} \phi \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}}+m_{i\eta} m_{j\eta} \frac{\partial^2 f_\mathrm{aniso}}{\partial m_{k\eta}\partial m_{\beta\eta}} \right]\right)\right\}
  \end{aligned}
\end{equation}

where $\phi$ is the shape function of the finite element method and if $k = \beta$, we have the on-diagonal contributions and if $k \neq \beta$ we have the off-diagonal contributions. Note that also, since $f_\mathrm{aniso}$ does not have cross-coupling between the sublattices (as is the case of other `Kernels` such as [AFMHomogeneousSublatticeExchange](source/kernels/AFMHomogeneousSublatticeExchange.md)), then further off-diagonal components ($\mathcal{J}_{m_{k\eta},m_{\beta\xi}}$) where $\eta \neq \xi$ are zeroed. As mentioned above, we have the `SC` suffix on this `Kernel` which means we are treating $\mathbf{A}$ as a MOOSE `Variable`. Therefore, we have to compute

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{m_{k\eta},A_{\beta}} &= \frac{\partial \mathcal{R}_{m_{k\eta}}}{\partial A_{\beta}} = \frac{\gamma}{\mu_0 M_s (1+\alpha^2)} \frac{\partial}{\partial A_{\beta}} \left\{\left(\psi_h, \epsilon_{ijk} m_{i\eta} \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} \right) + \alpha \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} m_{i\eta} m_{j\eta} \frac{\partial f_\mathrm{aniso}}{\partial m_{k\eta}} \right)\right\} \\
    &= \frac{\gamma}{\mu_0 M_s (1+\alpha^2)} \left\{\left(\psi_h, \epsilon_{ijk} m_{i\eta} \frac{\partial^2 f_\mathrm{aniso}}{\partial m_{k\eta} \partial A_\beta} \right) + \alpha \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} m_{i\eta} m_{j\eta} \frac{\partial^2 f_\mathrm{aniso}}{\partial m_{k\eta} \partial A_\beta} \right)\right\} \\
  \end{aligned}
\end{equation}

where only $\frac{\partial^2 f_\mathrm{aniso}}{\partial m_{k\eta}\partial A_\beta}$ is nonzero.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMSingleIonCubicSixthAnisotropySC object.

!syntax parameters /Kernels/AFMSingleIonCubicSixthAnisotropySC

!syntax inputs /Kernels/AFMSingleIonCubicSixthAnisotropySC

!syntax children /Kernels/AFMSingleIonCubicSixthAnisotropySC
