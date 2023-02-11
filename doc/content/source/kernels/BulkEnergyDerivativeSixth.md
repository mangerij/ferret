# BulkEnergyDerivativeSixth

!syntax description /Kernels/BulkEnergyDerivativeSixth

## Overview

Computes the residual and jacobian contributions for each component of the polarization vector $\mathbf{P}$. Note that if the structure as additional order parameters as is the case for perovskite $\mathrm{BiFeO}_3$ with its oxygen octahedral antiphase tilts, then this object can also be used for those components. The bulk free energy density, $f_\mathrm{bulk}$, up to sixth order is defined as,

\begin{equation}
  \begin{aligned}
    f_\mathrm{bulk} &= \alpha_{ij} P_i P_j + \alpha_{ijkl} P_i P_j P_k P_l + \alpha_{ijklmn} P_i P_j P_k P_l P_m P_n
  \end{aligned}
\end{equation}

Repeated indices are summed on for $i,j,k,l,... = x,y,z$. The relaxation time-evolution of the polarization $\mathbf{P}$ is given by,

\begin{equation}
  \begin{aligned}
    \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta f}{\delta \mathbf{P}},
  \end{aligned}
\end{equation}

so we need to compute the variational derivatives of $f$. Moving the RHS over, multiplying by a test function $\psi_h$ and integrating over the computational volume $\Omega$, we are left with the residual vector,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_\mathbf{P} = \left(\psi_h,\frac{\delta f}{\delta\mathbf{P}}\right).
  \end{aligned}
\end{equation}

The variational derivatives of $f_\mathbf{bulk}$ only involve derivatives with respect to each component so we can write explicitly the residual contribution as,

\begin{equation}
  \begin{aligned}
    \mathcal{R}_{P_\gamma} &= \left(\psi_h, 2 \alpha_{\gamma i} P_i + 4 \alpha_{\gamma ijk} P_i P_j P_k + ... \right), \\
  \end{aligned}
\end{equation}

and so on. To compute the on-diagonal jacobian contributions, we look for,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{P_\gamma, P_\gamma} = \frac{\partial \mathcal{R}_{P_\gamma}}{\partial P_\gamma} &= \left(\psi_h, 2 \phi \alpha_{\gamma \gamma} + 4 \phi \alpha_{\gamma \gamma ij} P_i P_j + ... \right), \\
  \end{aligned}
\end{equation}

with $\phi$ a shape function of the finite element method. The off-diagonal jacobian contributions provide the coupling to the other components of $\mathbf{P}$,

\begin{equation}
  \begin{aligned}
    \mathcal{J}_{P_\gamma, P_\beta} = \frac{\partial \mathcal{R}_{P_\gamma}}{\partial P_\beta} &= \left(\psi_h, 2 \phi \alpha_{\gamma \beta} + 4 \phi \alpha_{\gamma \beta ij} P_i P_j + ... \right). \\
  \end{aligned}
\end{equation}

If the material has a cubic parent phase (in the case of conventional perovskite ferroelectrics such as $\mathrm{PbTiO}_3, \mathrm{BaTiO}_3$, and $\mathrm{BiFeO}_3$), then most of these terms vanish and we are left with (in Voight notation) $\alpha_{1}, \alpha_{11}, \alpha_{12}, ..., \alpha_{123}$ nonzero. We should also note here in this overview that the sixth order expression for cubic phase materials has been 'hard-coded', but following our example other lower symmetry parent phases can be coded into FERRET.

## Example Input File Syntax

!! Describe and include an example of how to use the BulkEnergyDerivativeSixth object.

!syntax parameters /Kernels/BulkEnergyDerivativeSixth

!syntax inputs /Kernels/BulkEnergyDerivativeSixth

!syntax children /Kernels/BulkEnergyDerivativeSixth
