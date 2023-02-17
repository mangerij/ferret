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
    &= \frac{A_e}{4} \left\{\left(\nabla m_{1,x}\right)^2 + \left(\nabla m_{1,y}\right)^2 \left(\nabla m_{1,z}\right)^2 + \left(\nabla m_{2,x}\right)^2 + \left(\nabla m_{2,y}\right)^2 \left(\nabla m_{2,z}\right)^2\right\} + c.t.
  \end{aligned}
\end{equation}

TODO:

## Example Input File Syntax

!! Describe and include an example of how to use the AFMLocalSublatticeExchangeCartLL object.

!syntax parameters /Kernels/AFMLocalSublatticeExchangeCartLL

!syntax inputs /Kernels/AFMLocalSublatticeExchangeCartLL

!syntax children /Kernels/AFMLocalSublatticeExchangeCartLL
