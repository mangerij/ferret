# AFMLocalSublatticeExchangeCartLL

!syntax description /Kernels/AFMLocalSublatticeExchangeCartLL

## Overview

The long-range antiferromagnetic exchange free energy density can be written as

\begin{equation}
  \begin{aligned}
    f_\mathrm{\nabla L} = A_e \left[\left(\nabla L_x\right)^2 + \left(\nabla L_y\right)^2 \left(\nabla L_z\right)^2\right].
  \end{aligned}
\end{equation}

The quantity $\mathbf{L}$ is the antiferromagnetic Néel vector along with coupling constant $A_e$. Note that since $\mathbf{L} = \mathbf{m}_1 - \mathbf{m_2}/2$, we can write,

\begin{equation}
  \begin{aligned}
    f_\mathrm{\nabla L} &= A_e \left[\left(\nabla L_x\right)^2 + \left(\nabla L_y\right)^2 \left(\nabla L_z\right)^2\right] \\
    &= \frac{A_e}{4} \left\{\left(\nabla m_{1,x}\right)^2 + \left(\nabla m_{1,y}\right)^2 \left(\nabla m_{1,z}\right)^2 + \left(\nabla m_{2,x}\right)^2 + \left(\nabla m_{2,y}\right)^2 \left(\nabla m_{2,z}\right)^2\right\} \\ \nonumber
    &- 2 A_e \left(\frac{\partial m_{1,i}}{\partial x_j}\right)\left(\frac{\partial m_{2,i}}{\partial x_j}\right).
  \end{aligned}
\end{equation}

So we can separate this free energy density contribution into two distinct Kernels with the first being due to a standard expression in micromagnetics and the second being a cross-term that arises because of the coupling between the sublattices. The Landau-Lifshitz-Gilbert (LLG) equation for our two sublattice system is,

\begin{equation}\label{eqn:LLG_LLB}
  \begin{aligned}
    \frac{d \mathbf{m}_\eta}{dt} = -\frac{\gamma}{1+\alpha^2} \left(\mathbf{m}_\eta \times \mathbf{H}_\eta\right) - \frac{\gamma\alpha}{1+\alpha^2}\mathbf{m}_\eta\times\left(\mathbf{m}_\eta\times\mathbf{H}_\eta\right),
  \end{aligned}
\end{equation}

Multiplying by a test function $\psi_h$ suitable for finite element analysis, integrating over the volume, and neglecting the time depedent terms, we have in index notation,

\begin{equation}
  \begin{aligned}
    -\left(\psi_h, \frac{\gamma}{1+\alpha^2} \left(\mathbf{m}_\eta \times \mathbf{H}_\eta\right) \right) - \left(\psi_h, \frac{\gamma\alpha}{1+\alpha^2}\mathbf{m}_\eta\times\left(\mathbf{m}_\eta\times\mathbf{H}_\eta\right)\right) = 0
  \end{aligned}
\end{equation}


## Example Input File Syntax

!! Describe and include an example of how to use the AFMLocalSublatticeExchangeCartLL object.

!syntax parameters /Kernels/AFMLocalSublatticeExchangeCartLL

!syntax inputs /Kernels/AFMLocalSublatticeExchangeCartLL

!syntax children /Kernels/AFMLocalSublatticeExchangeCartLL
