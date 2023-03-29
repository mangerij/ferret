# AFMHomogeneousSublatticeExchange

!syntax description /Kernels/AFMHomogeneousSublatticeExchange

## Overview

Computes the residual and jacobian entries from the LLG equation due to a free energy density, $f_\mathrm{exch}$, corresponding to short-range sublattice exchange energy between two sublattice magnetizations $\{\mathbf{m}_1, \mathbf{m}_2\}$ in an antiferromagnet,

\begin{equation}
  \begin{aligned}
    f_\mathrm{exch} = 4 D_e \left(\mathbf{L}^2 - \mathbf{m}_2\right)
  \end{aligned}
\end{equation}

where $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$. Note that we do not use the conventional factor of 2 in these definitions but they can be evaluated after or during the simulation. Rewriting the above free energy, we have

\begin{equation}
  \begin{aligned}
    f_\mathrm{exch} = -4 D_e \left(m_{j1}\, m_{j2}\right)
  \end{aligned}
\end{equation}

The effective field due to this term can be calculated for each sublattice,

\begin{equation}
  \begin{aligned}
    H_{k1}^\mathrm{exch} &= \frac{4 D_e m_{k2}}{\mu_0 M_s}  \\
    H_{k2}^\mathrm{exch} &= \frac{4 D_e m_{k1}}{\mu_0 M_s}  \\
  \end{aligned}
\end{equation}

The effective field is evaluated at every time step of the two sublattice LLG-LLB equation,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{m}}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \mathbf{m}\times \mathbf{H} - \frac{\gamma \alpha}{1+\alpha^2} \mathbf{m} \times \mathbf{m}\times \mathbf{H},
  \end{aligned}
\end{equation}

where $\gamma$ is the electron gyromagnetic factor and $\alpha$ the phenomenological Gilbert damping parameter. Ignoring the time derivative, moving over the RHS to the LHS, and multiplying by a test function $\psi_h$, we have

\begin{equation}
  \begin{aligned}
    \frac{\gamma}{(1+\alpha^2)}\left(\psi_h, \mathbf{m}_\eta \times \mathbf{H}_\eta\right) + \frac{\gamma \alpha}{1+\alpha^2} \left(\psi_h, \mathbf{m}_\eta \times \mathbf{m}_\eta \times \mathbf{H}_\eta \right) = 0
  \end{aligned}
\end{equation}

Substituting the expression for $\mathbf{H}_\eta^\mathrm{exch}$ and writing the above in index notation gives the residual contribution for the $m_{k\eta}$ component,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{m_{k\eta}} = \frac{4 D_e \gamma}{\mu_0 M_s (1+\alpha^2)}\left(\psi_h, \epsilon_{ijk} m_{i\eta} m_{j\beta} \right) + \frac{4 D_e \gamma \alpha}{\mu_0 M_s \left(1+\alpha^2\right)} \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} m_{i\eta} m_{j\eta} m_{n\beta} \right),
  \end{aligned}
\end{equation}

where $\eta \neq \beta$ (i.e. $\eta = 1, \beta = 2$ or vice versa). To calculate the on- and off-diagonal jacobian entries, we need,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{m_{k\eta}, m_{\kappa\xi}} &= \frac{\partial \mathcal{R}_{m_{k\eta}}}{\partial m_{\kappa\xi}} = \frac{4 D_e \gamma}{\mu_0 M_s \left(1+\alpha^2\right)}\frac{\partial}{\partial m_{\kappa\xi} }\left\{ \left(\psi_h, \epsilon_{ijk} m_{i\eta} m_{j\beta} \right) + \alpha \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} m_{i\eta} m_{j\eta} m_{n\beta} \right)\right\} \\
    &= \frac{4 D_e \gamma}{\mu_0 M_s \left(1+\alpha^2\right)} \left\{ \left(\psi_h, \epsilon_{ijk} \left[ \delta_{\eta\xi} \delta_{i\kappa} m_{j\beta} + m_{i\eta} \delta_{\beta\xi} \delta_{j\kappa}\right] \right) + \alpha \left(\psi_h, \epsilon_{ijl}\epsilon_{lnk} \left[\delta_{i\kappa}\delta_{\eta\xi} m_{j\eta} m_{n\beta} + m_{i\eta} \delta_{j\kappa}\delta_{\eta\xi} m_{n\beta} + m_{i\eta} m_{j\eta} \delta_{n\kappa}\delta_{\eta\xi} \right] \right)\right\}\\
  \end{aligned}
\end{equation}

where $\xi = \eta$ and $\xi = \beta$ correspond to the on- and off-diagonal jacobian contributions.

## Example Input File Syntax

!! Describe and include an example of how to use the AFMHomogeneousSublatticeExchange object.

!syntax parameters /Kernels/AFMHomogeneousSublatticeExchange

!syntax inputs /Kernels/AFMHomogeneousSublatticeExchange

!syntax children /Kernels/AFMHomogeneousSublatticeExchange
