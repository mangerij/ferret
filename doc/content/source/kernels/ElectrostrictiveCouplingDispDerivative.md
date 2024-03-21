# ElectrostrictiveCouplingDispDerivative

!syntax description /Kernels/ElectrostrictiveCouplingDispDerivative

## Overview

We now consider the condition for mechanical equilibrium. The electric polarization, $\mathbf{P}$ generates a spontaneous contribution to the total strain,

\begin{equation}
  \begin{aligned}
    \varepsilon_{ij}^\mathrm{eig,P} = Q_{ijkl} P_k P_l,
  \end{aligned}
\end{equation}

where $Q_{ijkl}$ is the electrostrictive coefficient tensor. The governing equation of mechanical equilibrium is,

\begin{equation}
  \begin{aligned}
    \frac{\partial \sigma_{ij}}{\partial x_j} &= 0
    &= \frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + \varepsilon_{kl}^\mathrm{eig,P}\right)\right]
  \end{aligned}
\end{equation}

Multiplying both sides by the test function $\psi_h$ and integrating over the volumes,

\begin{equation}
  \begin{aligned}
    \left(\psi_h,\frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + \varepsilon_{kl}^\mathrm{eig,P}\right)\right]\right) &= 0 \\
    \left(\psi_h, C_{ijkl} \frac{\partial \varepsilon_{kl}}{\partial x_j}\right) + \left(\psi_h, C_{ijkl}\frac{\partial \varepsilon_{kl}^\mathrm{eig,P}}{\partial x_j}\right) &=0
  \end{aligned}
\end{equation}

Integrating by parts, we have,

\begin{equation}
  \begin{aligned}
    -\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}\right) + \left\langle\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}\right\rangle - \left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,P}\right) + \left\langle\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,P}\right\rangle  &=0.
  \end{aligned}
\end{equation}

The surface terms $\langle , \rangle$ vanish due to properties of the test function $\psi_h$ and we are left with the residual term for the elastic displacement component $u_i$,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{u_i} = -\underbrace{\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}\right)}_\mathrm{SolidMechanics} \,\,\,- \,\,\,\underbrace{\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl}^\mathrm{eig,P}\right)}_\mathrm{ElectrostrictiveCouplingDispDerivative}m
  \end{aligned}
\end{equation}

where the first term is handled by another `Kernel`. Using our expression for $\varepsilon_{kl}^\mathrm{eig,P}$, we find that this term does not explicitly depend on $\mathbf{u}$, therefore the on- and off-diagonal jacobian terms are

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{u_i, u_\beta} = \frac{\partial\mathcal{R}_{u_i}}{\partial u_\beta} = 0.
  \end{aligned}
\end{equation}

However, we still need to calculate,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{u_i, P_\beta} &= \frac{\partial\mathcal{R}_{u_i}}{\partial P_\beta} \\
    &= -\frac{\partial}{\partial A_\beta}\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} Q_{klmn} P_m P_n\right) \\
    &= - \left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} Q_{klmn} \phi \left[ \delta_{m\beta} P_n + P_m \delta_{n\beta}\right] \right)\\
    &= - \left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \phi \left[Q_{kl\beta n} P_n + Q_{klm\beta} P_m \right] \right)
  \end{aligned}
\end{equation}

where $\phi$ is the shape function of the finite element method.

## Example Input File Syntax

!! Describe and include an example of how to use the ElectrostrictiveCouplingDispDerivative object.

!syntax parameters /Kernels/ElectrostrictiveCouplingDispDerivative

!syntax inputs /Kernels/ElectrostrictiveCouplingDispDerivative

!syntax children /Kernels/ElectrostrictiveCouplingDispDerivative
