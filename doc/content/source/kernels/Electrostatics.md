# Electrostatics

!syntax description /Kernels/Electrostatics

## Overview

Laplace operator `Kernel` in the strong form is,

\begin{equation}
  \begin{aligned}
    \nabla A_\gamma \nabla \Phi_\mathrm{\gamma} = 0,
  \end{aligned}
\end{equation}

where $\gamma = \mathrm{H}$ or $\gamma = \mathrm{E}$ for magnetostatic or electrostatic calculations (or both). Repeated index summation is not implied in this `Kernel`. This means the material coefficient $A_\gamma = 1$ or $A_\gamma = \epsilon_b$ respectively. The parameter $\epsilon_b$ is the relative background dielectric permittivity in units of $\epsilon_0$. Multiplying by the test function $\psi_h$ and integrating both sides gives,

\begin{equation}
  \begin{aligned}
    \left(\psi_h , \nabla A_\gamma \nabla \Phi_\mathrm{\gamma}\right) = 0,
  \end{aligned}
\end{equation}

Integration by parts,

\begin{equation}
  \begin{aligned}
    -\left(A_\gamma \nabla \psi_h \cdot \nabla \Phi_\mathrm{\gamma}\right) + \left\langle A_\gamma \nabla \psi_h \cdot \nabla \Phi_\mathrm{\gamma}\right\rangle = 0,
  \end{aligned}
\end{equation}

Therefore, the residual contribution for $\Phi_\gamma$ is,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{\Phi_\gamma} = - \left(A_\gamma \nabla \psi_h \cdot \nabla \Phi_\mathrm{\gamma}\right),
  \end{aligned}
\end{equation}

with on-diagonal jacobian terms,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{\Phi_\gamma,\Phi_\gamma} = - \left(A_\gamma \nabla \psi_h \cdot \nabla \phi\right),
  \end{aligned}
\end{equation}

with $\phi$ the shape function of the finite element method.


## Example Input File Syntax

!! Describe and include an example of how to use the Electrostatics object.

!syntax parameters /Kernels/Electrostatics

!syntax inputs /Kernels/Electrostatics

!syntax children /Kernels/Electrostatics
