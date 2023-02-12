# WallEnergyDerivative

!alert construction title=Undocumented Class
The WallEnergyDerivative has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Kernels/WallEnergyDerivative

## Overview

This object computes the residual and jacobian contributions due to the variation of the free energy density associated with the Lifshitz invariants (gradients). We split the computation into two `Kernels`. The free energy density is,

\begin{equation}
  \begin{aligned}
  f_\mathrm{\nabla P} = \frac{G_{11}}{2}  \left( P_{x,x}^2 + P_{y,y}^2 + P_{z,z}^2 \right) +  G_{12}  \left(P_{x,x} P_{y,y} + P_{y,y} P_{z,z} + P_{x,x} P_{z,z} \right) + \frac{G_{44}}{2} \left[\left(P_{x,y} + P_{y,x} \right)^2+ \left(P_{y,z} + P_{z,y} \right)^2 + \left(P_{x,z} + P_{z,x}\right)^2\right]
  \end{aligned}
\end{equation}

 The governing equation of relaxation of the ferroelectric order is given by,

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


## Example Input File Syntax

!! Describe and include an example of how to use the WallEnergyDerivative object.

!syntax parameters /Kernels/WallEnergyDerivative

!syntax inputs /Kernels/WallEnergyDerivative

!syntax children /Kernels/WallEnergyDerivative
