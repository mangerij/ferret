# AFMInteractionCartLLHConst

!alert construction title=Undocumented Class
The AFMInteractionCartLLHConst has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Kernels/AFMInteractionCartLLHConst

## Overview

Consider the Zeeman interaction energy density,

\begin{equation}
  \begin{aligned}
    f_\mathrm{Zeeman} = - M_s \mathbf{m} \cdot \mathbf{H}_\mathrm{appl}
  \end{aligned}
\end{equation}

where, in an antiferromagnet (AFM), the total net magnetization is $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$. Therefore, the effective field due to this term is,

\begin{equation}
  \begin{aligned}
    \mathbf{H}_{\mathrm{appl},\eta} = - \frac{1}{\mu_0 M_s} \frac{\delta f_\mathrm{Zeeman}}{\delta \mathbf{m}_\eta}
  \end{aligned}
\end{equation}

with $\eta = 1,2$ being the index of the sublattice. We aim to solve the normalized two-sublattice LLG equation,

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

To be finshed..

## Example Input File Syntax

!! Describe and include an example of how to use the AFMInteractionCartLLHConst object.

!syntax parameters /Kernels/AFMInteractionCartLLHConst

!syntax inputs /Kernels/AFMInteractionCartLLHConst

!syntax children /Kernels/AFMInteractionCartLLHConst
