!include tutorials/tutorial_header.md

# Tutorial 3: Ferromagnetic ringdown

This tutorial (and others) covers the basic usage of the micromagnetics implemented in FERRET. This simple problem is for solving the time-dependent LLG equation with exchange and demagnetizing field phenomena. We consider a magnetic body with a geometry $(20\times 20\times 3)$ which we define with the Mesh block.

!listing tutorial/ringdown.i
         block=Mesh
         link=False
         language=python

In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through Materials, Kernels, or other MOOSE objects. More on this later in this tutorial. Note that the total computational domain is $(50\times 50\times 20)$. The remaining area is defined as a vacuum block which acts as a numerical resource in order to solve for the demagnetizing field far from the magnetic body. In principle, we should take the $r\to \infty$ limit which is impossible in a numerical simulation so this area is used must be large enough such that the field from the magnetic body goes to zero.

We use the SubdomainBoundingBoxGenerator to split the mesh into two blocks. In this problem,

!listing tutorial/ringdown.i
         block=Variables
         link=False
         language=python

shows that we are only solving for the variables $(m_x, m_y, m_y$, and $\Phi_\mathrm{H})$. The initial conditions (ICs) are seeded via spherical coordinates using RandomConstrainedVectorFieldIC to ensure that the ICs satisfy the normalization condition of the LLG equation, $|\mathbf{m}| = 1$.

The micromagnetic problem solves the LLG equation,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{M}}{\partial t} = -\gamma \mathbf{M}\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M} \times \mathbf{M}\times \mathbf{H}
  \end{aligned}
\end{equation}

which is partitioned into a few Kernel objects for computation. We separate contributions component-wise and also by effective fields due to different physics in the problem (in this case exchange stiffness and demagnetization phenomena).

These Kernels are listed

!listing tutorial/ringdown.i
         block=AuxKernels
         link=False
         language=python

!listing tutorial/ringdown.i
         block=Kernels
         link=False
         language=python

where

The Materials block,

!listing tutorial/ringdown.i
         block=Materials
         link=False
         language=python

shows the coefficients we will use. Some comments here are also warranted.
