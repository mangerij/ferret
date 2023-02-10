# ConversePiezoelectricStrain

!alert construction title=Undocumented Class
The ConversePiezoelectricStrain has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Kernels/ConversePiezoelectricStrain

## Overview

In the linear limit, the constitutive governing equation for piezoelectricity is

\begin{equation}
  \begin{aligned}
    \frac{\partial}{\partial x_j} \left(C_{ijkl} \varepsilon_{kl} + d_{klm} E_m\right) = 0,
  \end{aligned}
\end{equation}

where $C_{ijkl}, \varepsilon_{kl},$ and $E_m$ are the components of the elastic stiffness tensor of rank four, elastic strain tensor of rank two, and the electric field vector components. The direct piezoeletric coefficient $d_{klm}$ is of rank three. This is equivalently the equation for mechanical equilibrium with coupling to the electric field,  $E_m = - \nabla \Phi_\mathrm{E}$. Note that the strain tensor is defined in the usual way, $\varepsilon_{kl} = \frac{1}{2} \left( \frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k} \right)$. Einstein summation is implied. Multiplying both sides by a test function $\psi_h$ and integrating over the volume yields

\begin{equation}
  \begin{aligned}
    \left(\psi_h, C_{ijkl} \frac{\partial \varepsilon_{kl} }{\partial x_j} \right) + \left(\psi_h, d_{klm} \frac{\partial E_m}{\partial x_j} \right) = 0.
  \end{aligned}
\end{equation}

Using integration by parts and substituting $E_m = - \partial \Phi_\mathrm{E} / \partial x_m$ gives,

\begin{equation}
  \begin{aligned}
    \underbrace{-\left(\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl} \right) + \left\langle\frac{\partial \psi_h}{\partial x_j}, C_{ijkl} \varepsilon_{kl} \right\rangle}_{\mathrm{TensorMechanics}} \,\,\,\,\, + \,\,\,\,\, \underbrace{\left( \frac{\partial \psi_h}{\partial x_j}, d_{klm} \frac{\partial \Phi_\mathrm{E}}{\partial x_m}\right) - \left\langle \frac{\partial \psi_h}{\partial x_j}, d_{klm} \frac{\partial \Phi_\mathrm{E}}{\partial x_m}\right\rangle}_{\mathrm{ConversePiezoelectricStrain}} + s.t. = 0.
  \end{aligned}
\end{equation}

Due to the properties of the test function, the surface terms $\langle . \rangle$ vanish and we are left with the residual contribution due to the converse piezoelectric strain. Note that the `TensorMechanics` `Action` system automatically sets up and handles the first term in the above equation. Therefore, one needs `TensorMechanics` active in the `Kernels` block of their input files.

## Example Input File Syntax

!! Describe and include an example of how to use the ConversePiezoelectricStrain object.

!syntax parameters /Kernels/ConversePiezoelectricStrain

!syntax inputs /Kernels/ConversePiezoelectricStrain

!syntax children /Kernels/ConversePiezoelectricStrain
