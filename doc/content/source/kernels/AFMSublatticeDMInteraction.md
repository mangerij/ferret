# AFMSublatticeDMInteraction

!syntax description /Kernels/AFMSublatticeDMInteraction

## Overview

Consider the free energy density contribution due to the Dzyaloshinskii–Moriya interaction (DMI),

\begin{equation}
  \begin{aligned}
    f_\mathrm{DMI} = 4 D_0 \mathbf{A} \cdot \left( \mathbf{L} \times \mathbf{m}\right)
  \end{aligned}
\end{equation}

where $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$. The effective field due to this term is,

\begin{equation}
  \begin{aligned}
    H_{k\eta}^\mathrm{DMI} &= - \frac{1}{\mu_0 M_s} \frac{\delta f_\mathrm{DMI}}{\delta m_{k\eta}}.\\
  \end{aligned}
\end{equation}

We can write

\begin{equation}
  \begin{aligned}
    f_\mathrm{DMI} = 4 D_0 A_k \epsilon_{ijk} \left(m_{i1} m_{j1} - m_{i2} m_{i1} + m_{i1} m_{j2} - m_{i2} m_{j2}\right)
  \end{aligned}
\end{equation}

which implies,

\begin{equation}
  \begin{aligned}
    H_{k1}^\mathrm{DMI} &= - \frac{1}{\mu_0 M_s} \left\{\delta_{ik} m_{j1} + m_{i1} \delta_{jk} - m_{i2} \delta){ik} + \delta_{ik} m_{j2}\right\}.\\
    H_{k2}^\mathrm{DMI} &= - \frac{1}{\mu_0 M_s} \left\{- \delta_{ik} m_{i1} + m_{i1} \delta_{jk} - \delta_{ik} m_{j2} - m_{i2} \delta_{jk}\right\}
  \end{aligned}
\end{equation}

explicitly for the two sublattices $\eta = 1,2$. The effective field is evaluated at every time step of the two sublattice LLG-LLB equation,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{m}_\eta}{\partial t} = -\frac{\gamma}{(1+\alpha^2)} \mathbf{m}_\eta\times \mathbf{H}_\eta - \frac{\gamma \alpha}{1+\alpha^2} \mathbf{m}_\eta \times \mathbf{m}_\eta\times \mathbf{H}_\eta,
  \end{aligned}
\end{equation}

Substituting our expressions for $H_{k\eta}^\mathrm{DMI}$, we have, 


## Example Input File Syntax

!! Describe and include an example of how to use the AFMSublatticeDMInteraction object.

!syntax parameters /Kernels/AFMSublatticeDMInteraction

!syntax inputs /Kernels/AFMSublatticeDMInteraction

!syntax children /Kernels/AFMSublatticeDMInteraction
