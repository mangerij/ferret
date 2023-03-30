# AFMEasyPlaneAnisotropy

!syntax description /Kernels/AFMEasyPlaneAnisotropy

## Overview

Computes the residual and jacobian entries from the LLG equation due to a free energy density, $f_\mathrm{easy-plane}$, corresponding to easy-plane magnetocrystalline anisotropy,

\begin{equation}
  \begin{aligned}
    f_\mathrm{easy-plane} = \sum\limits_{\eta = 1,2} K_1 \left(\textbf{m}_\eta \cdot \hat{\mathbf{P}}\right)^2
  \end{aligned}
\end{equation}

where $\hat{\mathbf{P}}$ is the director of the electric polarization field and sublattice indices $\eta = 1,2$. The effective field due to this term can be calculated, in index notation,

\begin{equation}
  \begin{aligned}
    H_{k\eta}^\mathrm{easy-plane} &= - \frac{1}{\mu_0 M_s} \frac{\delta f_\mathrm{easy-plane}}{\delta m_{k\eta}}\\
    &= - \frac{2}{\mu_0 M_s} \hat{P}_k (m_{j\eta} \hat{P}_j)^2 \\
  \end{aligned}
\end{equation}

The effective field is evaluated at every time step of the two sublattice LLG-LLB equation,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{m}_\eta}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \mathbf{m}_\eta\times \mathbf{H}_\eta - \frac{\gamma \alpha}{1+\alpha^2} \mathbf{m}_\eta \times \mathbf{m}_\eta\times \mathbf{H}_\eta,
  \end{aligned}
\end{equation}

where $\gamma$ is the electron gyromagnetic factor and $\alpha$ the phenomenological Gilbert damping parameter. Ignoring the time derivative, moving over the RHS to the LHS, and multiplying by a test function $\psi_h$, we have

\begin{equation}
  \begin{aligned}
    \frac{\gamma}{(1+\alpha^2)}\left(\psi_h, \mathbf{m}_\eta \times \mathbf{H}_\eta\right) + \frac{\gamma \alpha}{1+\alpha^2} \left(\psi_h, \mathbf{m}_\eta \times \mathbf{m}_\eta \times \mathbf{H}_\eta \right) = 0
  \end{aligned}
\end{equation}

In index notation, we have, with $\eta$ a free index,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_{n\eta}} = \frac{\gamma}{(1+\alpha^2)}\left(\psi_h, \epsilon_{ijn} m_{i\eta} H_{j\eta}\right) + \frac{\gamma \alpha}{1+\alpha^2} \left(\psi_h, \epsilon_{ijk} \epsilon_{kln} m_{i\eta} m_{j\eta}  H_{l\eta} \right)
  \end{aligned}
\end{equation}

Substituting our definitions of $H_{n\eta}$ we have the residual contribution for the components of $\mathbf{m}_\eta$,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_{n\eta}} = - \frac{2 \gamma}{\mu_0 M_s (1+\alpha^2)} \left(\psi_h, \epsilon_{ijn} m_{i\eta} \hat{P}_j (m_{q\eta} \hat{P}_q)^2 \right) - \frac{2 \gamma \alpha}{\mu_0 M_s \left(1+\alpha^2\right)} \left(\psi_h, \epsilon_{ijk} \epsilon_{kln} m_{i\eta} m_{j\eta}  \hat{P}_l (m_{q\eta} \hat{P}_q)^2 \right).
  \end{aligned}
\end{equation}

The jacobian contributions are calculated by,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{m_{n\eta},m_{\beta\eta}} &= \frac{\partial \mathcal{R}_{m_{n\eta}}}{\partial m_{\beta\eta}} \\
    &= \frac{\partial}{\partial m_{\beta\eta}} \left\{- \frac{2 \gamma}{\mu_0 M_s (1+\alpha^2)} \left(\psi_h, \epsilon_{ijn} m_{i\eta} \hat{P}_j (m_{q\eta} \hat{P}_q)^2 \right) - \frac{2 \gamma \alpha}{\mu_0 M_s \left(1+\alpha^2\right)} \left(\psi_h, \epsilon_{ijk} \epsilon_{kln} m_{i\eta} m_{j\eta}  \hat{P}_l (m_{q\eta} \hat{P}_q)^2 \right)\right\}. \\
    &= - \frac{2 \gamma}{\mu_0 M_s (1+\alpha^2)} \left( \psi_h, \epsilon_{ijn} \hat{P}_j \hat{P}_q m_{q\eta} \phi \left\{\delta_{i\beta} \left(m_{r\eta} \hat{P}_r\right) + 2 m_{i\eta} \hat{P}_\beta \right\} \right) - \frac{2 \gamma \alpha}{\mu_0 M_s (1+\alpha^2)} \left( \psi_h, \epsilon_{ijk} \epsilon_{kln} \hat{P}_l  m_{q\eta} \hat{P}_q \phi \left\{ m_{r\eta} \hat{P}_r \left(\delta_{i\beta} + \delta_{j\beta}\right) + 2 m_{i\eta} m_{j\eta} \right\} \right)
  \end{aligned}
\end{equation}

where we use $\partial m_{j\eta} / \partial m_{k\eta} = \delta_{jk} \phi$ with $\phi$ a shape function of the finite element method.


## Example Input File Syntax

!! Describe and include an example of how to use the AFMEasyPlaneAnisotropy object.

!syntax parameters /Kernels/AFMEasyPlaneAnisotropy

!syntax inputs /Kernels/AFMEasyPlaneAnisotropy

!syntax children /Kernels/AFMEasyPlaneAnisotropy
