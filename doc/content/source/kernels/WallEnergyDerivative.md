# WallEnergyDerivative

!alert construction title=Undocumented Class
The WallEnergyDerivative has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Kernels/WallEnergyDerivative

## Overview

This object computes the residual and jacobian contributions due to the variation of the free energy density associated with the Lifshitz invariants (gradients). We split the computation into two `Kernels`. The gradient free energy density for a parent phase cubic perovskite is (in Voight notation),

\begin{equation}
  \begin{aligned}
  f_\mathrm{\nabla P} = \frac{G_{11}}{2}  \left( P_{x,x}^2 + P_{y,y}^2 + P_{z,z}^2 \right) +  G_{12}  \left(P_{x,x} P_{y,y} + P_{y,y} P_{z,z} + P_{x,x} P_{z,z} \right) + \frac{G_{44}}{2} \left[\left(P_{x,y} + P_{y,x} \right)^2+ \left(P_{y,z} + P_{z,y} \right)^2 + \left(P_{x,z} + P_{z,x}\right)^2\right]
  \end{aligned}
\end{equation}

which can be written in general as,

\begin{equation}
  \begin{aligned}
  f_\mathrm{\nabla P} = G_{ijkl}\frac{\partial P_k}{\partial x_l}.
  \end{aligned}
\end{equation}

The governing time-dependent Landau-Ginzburg-Devonshire (TDLGD) equation of relaxation of the ferroelectric order is given by,

\begin{equation}
 \begin{aligned}
   \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta f}{\delta \mathbf{P}},
 \end{aligned}
\end{equation}

which means we need to compute $\frac{\delta f}{\delta \mathbf{P}}$. Multiplying by the test function $\psi_h$ and integrating over the volume, we have after moving $\Gamma_P$ over to the other side,

\begin{equation}
 \begin{aligned}
   \left(\psi_h,\frac{1}{\Gamma_P}\frac{\partial \mathbf{P}}{\partial t}\right) + \left(\psi_h, \frac{\delta f}{\delta \mathbf{P}} \right) = 0
 \end{aligned}
\end{equation}

The variational derivative of the free energy is as follows,

\begin{equation}
 \begin{aligned}
   \frac{\delta f}{\delta \mathbf{P}} = \frac{\partial f}{\partial \mathbf{P}} - \frac{\partial}{\partial \mathbf{r}} \cdot \frac{\partial f}{\left(\frac{\partial \mathbf{P}}{\partial \mathbf{r}}\right)}
 \end{aligned}
\end{equation}

Focusing on $f_\mathrm{\nabla P}$, neglecting the time derivative, and writing in index notation we have

\begin{equation}
 \begin{aligned}
   \left(\psi_h,\frac{\partial f_\mathrm{\nabla P}}{\partial P_\beta}\right) - \frac{\partial}{\partial x_\eta} \left(\psi_h, \frac{\partial f_\mathrm{\nabla P}}{\partial \left(\frac{\partial P_\beta}{\partial x_\eta}\right)}\right) &= 0 \\
   \left(\psi_h,\frac{\partial G_{ijkl} \frac{\partial P_k}{\partial x_l}  }{\partial P_\beta}\right) - \frac{\partial}{\partial x_\eta} \left(\psi_h, \frac{\partial  G_{ijkl} \frac{\partial P_k}{\partial x_l}}{\partial \left(\frac{\partial P_\beta}{\partial x_\eta}\right)}\right) &= 0.
 \end{aligned}
\end{equation}

In the first term, we use the relation

\begin{equation}
 \begin{aligned}
   \frac{\partial \left(G_{ijkl} \frac{\partial P_k}{\partial x_l}\right)}{\partial P_\beta} &= \delta_{k\beta} 4 G_{ijkl} \frac{\partial \phi}{\partial x_l} \\
   &= 4 G_{ij\beta l} \frac{\partial \phi}{\partial x_l}
 \end{aligned}
\end{equation}

with $\phi$ the shape function of the finite element method. The second term of $\left(\psi_h, delta f_\mathrm{\nabla P} / \delta \mathbf{P}\right)$ can be written as

\begin{equation}
 \begin{aligned}
   - \frac{\partial}{\partial x_\eta} \left(\psi_h, \frac{\partial  \left(G_{ijkl} \frac{\partial P_k}{\partial x_l}\right)}{\partial \left(\frac{\partial P_\beta}{\partial x_\eta}\right)}\right) &= \left(\frac{\partial \psi_h}{\partial x_\eta}, \frac{\partial  \left(G_{ijkl} \frac{\partial P_k}{\partial x_l}\right)}{\partial \left(\frac{\partial P_\beta}{\partial x_\eta}\right)}\right) - \left\langle \frac{\psi_h}{\partial x_\eta}, \frac{\partial  \left(G_{ijkl} \frac{\partial P_k}{\partial x_l}\right)}{\partial \left(\frac{\partial P_\beta}{\partial x_\eta}\right)}\right\rangle\\
 \end{aligned}
\end{equation}

via integration by parts. The surface terms vanish and we are left with the residual contribution for $P_\beta$ due to the Lifshitz invariants,

\begin{equation}
 \begin{aligned}
  \mathcal{R}_{P_\beta} = \underbrace{\left(\frac{\partial \psi_h}{\partial x_\eta}, \frac{\partial  \left(G_{ijkl} \frac{\partial P_k}{\partial x_l}\right)}{\partial \left(\frac{\partial P_\beta}{\partial x_\eta}\right)}\right)}_\mathrm{WallEnergyDerivative} \,\,\,\,+\,\,\,\, \underbrace{\left(\psi_h, 4 G_{ij\beta l} \frac{\partial \phi}{\partial x_l}\right)}_\mathrm{Wall2EnergyDerivative} \\
 \end{aligned}
\end{equation}

where the first term is computed in `WallEnergyDerivative` and the second is calculated in `Wall2EnergyDerivative`. We hand code $G_{ijkl}$ into the problem for cubic parent phase perovskites such as $\mathrm{BaTiO}_3$. This term can be reworked in general for lower symmetry ferroelectric materials if needed. In general, due to the fact that the second term scales on the order of $\psi_h \partial \phi / \partial x_k$, its influence to the evolution of the TDLGD problem is quite weak. Now we compute the jacobian contributions,

\begin{equation}
 \begin{aligned}
  \mathcal{J}_{P_\beta, P_xi} = \frac{\partial \mathcal{R}_{P_\beta}}{\partial P_\xi} =  \\
 \end{aligned}
\end{equation}

## Example Input File Syntax

!! Describe and include an example of how to use the WallEnergyDerivative object.

!syntax parameters /Kernels/WallEnergyDerivative

!syntax inputs /Kernels/WallEnergyDerivative

!syntax children /Kernels/WallEnergyDerivative
