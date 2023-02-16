# Constitutive theory of piezoelectric materials

Within FERRET/MOOSE we have implemented the governing equations of piezoelectrics which allows the user to simulate strongly coupled electromechanical phenomena. For a linear piezoelectric material, the stress-divergence equation for mechanical equilibrium reads

\begin{equation}
  \begin{aligned}
    \frac{\partial}{\partial x_j} \left[C_{ijkl} \left()\varepsilon_{kl} + d_{klm} E_m\right)\right] = 0,
  \end{aligned}
\end{equation}

where $C_{ijkl}, \varepsilon_{kl},$ and $E_m$ are the components of the elastic stiffness tensor, elastic strain tensor, and electric field. The direct piezoeletric coefficient $d_{klm}$ is of rank three. Note that the strain tensor is defined in the usual way, $\varepsilon_{kl} = \frac{1}{2} \left( \frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k} \right).$

Additionally, the Poisson equation also includes contributions from the (converse piezo) strain-charge,

\begin{equation}
  \begin{aligned}
    \frac{\partial^2 \Phi_\mathrm{E} }{\partial x_r \partial x_r} + \frac{\partial d_{jkl} \sigma_{kl}}{\partial x_j} = 0 ,
  \end{aligned}
\end{equation}

 with $\sigma_{kl}$ being the components of the elastic stress tensor. With sufficient choice of materials parameters, Equations (1) and (2) can be solved self-consistently under arbitrary mechanical loads or applied electric fields to yield the static configuration of the electrostatic potential and elastic displacements $\mathbf{u}.$

Equations (1) and (2) can be cast dynamically, to simulate piezeoelectric actuation in real time. This is done by setting the LHS of the first equation to be equal to $\partial u_i / \partial t$. The time scale can be set by a prefactor.


Another possible use of this implementation is to rotate all of the tensorial coefficients grain-wise leading to calculations of a polycrystalline piezoelectric. Within our block-restricted polycrystal approach, this allows for the piezoelectric properties to be evaluated in a computational box with a \textit{real} grain structure by rotating the tensors via,

\begin{equation}
  \begin{aligned}
    \tilde{A}_{ij} &= R_{i\alpha}R_{j\beta} A_{\alpha \beta} \\
    \tilde{A}_{ijk} &= R_{i\alpha}R_{j\beta}R_{k\eta} A_{\alpha \beta \eta} \\
  \end{aligned}
\end{equation}

where $A_{ij}$ and $A_{ijk}$ are second and third rank tensors and $R_{ij}$ is a rotation operator using the internal `RotationTensor` operation in MOOSE utils. The rotation operator $R_{ij}$ accepts Euler angles in the standard Bunge sequence ($\mathbf{ZXZ}$).
