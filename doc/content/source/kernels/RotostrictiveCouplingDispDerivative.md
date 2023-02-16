# RotostrictiveCouplingDispDerivative

!syntax description /Kernels/RotostrictiveCouplingDispDerivative

## Overview

We now consider the condition for mechanical equilibrium. In addition to the electric polarization, $\mathbf{P}$, there are ferroelectrics where secondary structural order parameters condense below the phase transition temperature. In the case of $\mathrm{BiFeO}_3$, the oxygen octahedral tilt cages can generate a contribution to the total spontaneous strain,

\begin{equation}
  \begin{aligned}
    \varepsilon_{ij}^\mathrm{eig,A} = R_{ijkl} A_k A_l,
  \end{aligned}
\end{equation}

where $R_{ijkl}$ is the rotostrictive coefficient tensor and $\mathbf{A}$ is the order parameter associated with the antiphase tilts. The governing equation of mechanical equilibrium is,

\begin{equation}
  \begin{aligned}
    \frac{\partial \sigma_{ij}}{\partial x_j} &= 0
    &= \frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + \varepsilon_{kl}^\mathrm{eig,P} + \varepsilon_{kl}^\mathrm{eig,A}\right)\right]
  \end{aligned}
\end{equation}

Multiplying both sides by the test function $\psi_h$ and integrating over the volumes,

\begin{equation}
  \begin{aligned}
    \left(\psi_h,\frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + \varepsilon_{kl}^\mathrm{eig,P} + \varepsilon_{kl}^\mathrm{eig,A}\right)\right] \right) &= 0 \\
    \left(\psi_h, C_{ijkl} \frac{\partial \varepsilon_{kl}}{\partial x_j}\right) + \left(\psi_h, C_{ijkl}\frac{\partial \varepsilon_{kl}^\mathrm{eig,P}}{\partial x_j}\right) + \left(\psi_h, C_{ijkl}\frac{\partial \varepsilon_{kl}^\mathrm{eig,A}}{\partial x_j}\right) &=0
  \end{aligned}
\end{equation}

Integrating by parts, we have,

\begin{equation}
  \begin{aligned}
    -\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}\right) + \left\langle\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}\right\rangle - \left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,P}\right) + \left\langle\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,P}\right\rangle -\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,A}\right) + \left\langle\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,A}\right\rangle &=0.
  \end{aligned}
\end{equation}

The surface terms $\langle , \rangle$ vanish due to properties of the test function $\psi_h$ and we are left with the residual term for the elastic displacement component $u_i$,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{u_i} = -\underbrace{\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}\right)}_\mathrm{TensorMechanics} \,\,\,- \,\,\,\underbrace{\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,P}\right)}_\mathrm{ElectrostrictiveCouplingDispDerivative} \,\,\, - \,\,\,\underbrace{\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,A}\right)}_\mathrm{RotostrictiveCouplingDispDerivative}.
  \end{aligned}
\end{equation}

where the first two terms are handled by other `Kernels`. Using our expression for $\varepsilon_{kl}^\mathrm{eig,A}$, we find that this term does not explicitly depend on $\mathbf{u}$, therefore the on- and off-diagonal jacobian terms are

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{u_i, u_\beta} = \frac{\partial\mathcal{R}_{u_i}}{\partial u_\beta} = 0.
  \end{aligned}
\end{equation}

However, we still need to calculate,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{u_i, A_\beta} &= \frac{\partial\mathcal{R}_{u_i}}{\partial A_\beta} \\
    &= -\frac{\partial}{\partial A_\beta}\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} R_{klmn} A_m A_n\right) \\
    &= - \left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} R_{klmn} \phi \left[ \delta_{m\beta} A_n + A_m \delta_{n\beta}\right] \right)\\
    &= - \left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \phi \left[R_{kl\beta n} A_n + R_{klm\beta} A_m \right] \right)
  \end{aligned}
\end{equation}

where $\phi$ is the shape function of the finite element method.



## Example Input File Syntax

!! Describe and include an example of how to use the RotostrictiveCouplingDispDerivative object.

!syntax parameters /Kernels/RotostrictiveCouplingDispDerivative

!syntax inputs /Kernels/RotostrictiveCouplingDispDerivative

!syntax children /Kernels/RotostrictiveCouplingDispDerivative
