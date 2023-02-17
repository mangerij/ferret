# SurfaceMechanicsBC

!syntax description /BCs/SurfaceMechanicsBC

## Overview

An `IntegratedBC` object which computes a residual contribution due to a surface tension $\tau$. The surface elastic energy density is defined as,

\begin{equation}
  \begin{aligned}
    F_\mathrm{surface} = \int\limits_S d^2 \mathbf{r} \sigma_{\alpha \beta}^s \varepsilon_{\alpha \beta}
  \end{aligned}
\end{equation}

where $\alpha,\beta = 1,2$ denote a local orthonormal basis in the tangent plane of an arbitrary surface $S$ with a surface normal $\mathbf{n}(\mathbf{r})$. The surface stress is related to the surface elastic stiffness tensor via,

\begin{equation}
  \begin{aligned}
    \sigma_{\alpha\beta} &= \tau_{\alpha\beta} + C_{\alpha \beta \delta \gamma}^s \varepsilon_{\delta \gamma}^s \\
    &= \tau \delta_{\alpha\beta} + C_{\alpha \beta \delta \gamma}^s \varepsilon_{\delta \gamma}^s
  \end{aligned}
\end{equation}

to first approximation ($\tau$) being an average of the residual surface stress components and $\delta_{\alpha\beta}$ the Kronecker symbol. A projection operator can be constructed such that

\begin{equation}
  \begin{aligned}
    \hat{\mathbf{P}} &= \mathcal{1} - \mathbf{n} \outerp \mathbf{n}
  \end{aligned}
\end{equation}

Therefore, $\sigma_{ij}^s = \hat{P}_{i\alpha} \hat{P}_{j\beta} \sigma_{\alpha \beta}$ and $\varepsilon_{ij}^s = \hat{P}_{i\alpha} \hat{P}_{j\beta} \varepsilon_{\alpha \beta}$.
Now we need to consider the stress divergence equation (condition of mechanical equilibrium of the surface) which is

\begin{equation}
  \begin{aligned}
    \frac{\partial \sigma_{\alpha \beta}^2}{\partial x_\beta} &= 0\\
    &= \frac{\partial}{\partial x_\beta} \left(\tau \delta_{\alpha\beta} + C_{\alpha \beta \delta \gamma}^s \varepsilon_{\delta \gamma}^s\right)\\
    &= \frac{\partial}{\partial x_\beta} \left(\tau \delta_{\alpha\beta} + C_{\alpha \beta \delta \gamma}^s \hat{P}_{i\delta} \hat{P}_{j\gamma} \varepsilon_{\delta\gamma}\right)
  \end{aligned}
\end{equation}

TODO:


## Example Input File Syntax

!! Describe and include an example of how to use the SurfaceMechanicsBC object.

!syntax parameters /BCs/SurfaceMechanicsBC

!syntax inputs /BCs/SurfaceMechanicsBC

!syntax children /BCs/SurfaceMechanicsBC
