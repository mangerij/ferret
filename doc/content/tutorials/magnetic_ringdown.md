!include tutorials/tutorial_header.md

# Tutorial 4: Ferromagnetic ringdown

This tutorial (and others) covers the basic usage of the micromagnetics implemented in FERRET. This simple problem is for solving the time-dependent LLG equation with exchange stiffness and demagnetizing field phenomena in a permalloy magnet. The total energy density of the system is written as

\begin{equation}
  \begin{aligned}
    f_\mathrm{total} = f_\mathrm{exch} + f_\mathrm{magnetostatic}
  \end{aligned}
\end{equation}

where $f_\mathrm{exch} = A_e \mathbf{m}\cdot \nabla^2 \mathrm{m}$ and $f_\mathrm{magnetostatic} = -\mathbf{M}\cdot\mathbf{H} = M_s \mathbf{m}\cdot \Phi_\mathrm{H}$. The coefficient $A_e$ is the exchange stiffness parameter, $M_s$ the saturation magnetization density, and $\Phi_\mathrm{H}$ the magnetostatic potential.

We consider a magnetic body with a geometry $(20\times 20\times 3)$ which we define with the `Mesh` block.

!listing tutorial/ringdown.i
         block=Mesh
         link=False
         language=python

In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through Materials, Kernels, or other MOOSE objects. More on this later in this tutorial. Note that the total computational domain is $(50\times 50\times 20)$. The remaining area is defined as a vacuum block which acts as a numerical resource in order to solve for the demagnetizing field far from the magnetic body. In principle, we should take the $r\to \infty$ limit which is impossible in a numerical simulation so the volume that is used must be large enough such that the field from the magnetic body naturally goes to zero. We suggest that if the vacuum block approach is used, one should likely test a number of distances to ensure the expected $1/r^3$ dependence of the demagnetizing field. We use the MOOSE object `SubdomainBoundingBoxGenerator` to split the mesh into two blocks. In this problem,

!listing tutorial/ringdown.i
         block=Variables
         link=False
         language=python

shows that we are only solving for the variables $(m_x, m_y, m_y$, and $\Phi_\mathrm{H})$. The initial conditions (ICs) are seeded via spherical coordinates using [`RandomConstrainedVectorFieldIC`](source/ics/RandomConstrainedVectorFieldIC.md) to ensure that the ICs satisfy the normalization condition of the LLG equation, $|\mathbf{m}| = 1$.

The micromagnetic problem solves the LLG-LLB equation,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{M}}{\partial t} = -\gamma \mathbf{M}\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M} \times \mathbf{M}\times \mathbf{H} - \frac{\gamma \alpha_\mathrm{LLB}}{1+\alpha^2} \left(m^2 - 1\right) \mathbf{m}
  \end{aligned}
\end{equation}

with $\mathbf{H} = - (1 / \mu_0) \delta F /\delta \mathbf{M}$. The LLG-LLB equation is partitioned into a few Kernel objects for computation. We separate contributions component-wise and also by specific effective fields (in this case exchange stiffness and demagnetizing fields). These Kernels are listed

!listing tutorial/ringdown.i
         block=Kernels
         link=False
         language=python

where [`MasterExchangeCartLLG`](source/kernels/MasterExchangeCartLLG.md) handles the terms involving the exchange stiffness energy density, [`MasterInteractionCartLLG`](source/kernels/MasterInteractionCartLLG.md) involves the interaction with the demagnetizing or applied magnetic fields, [`MasterLongitudinalLLB`](source/kernels/MasterLongitudinalLLB.md) is the longitudinal restoring term from the LLB approximation. Finally, we have  [`MagHStrongCart`](source/kernels/MagHStrongCart.md) and [`Electrostatics`](source/kernels/Electrostatics.md) (the Laplace equation LHS $\nabla^2 \Phi_\mathrm{H}$) for the magnetostatic Poisson equation solved at every time step. The user can click these hyperlinks to see the weak-form algebra necessary to construct the `Kernel` objects. The `Materials` block,

!listing tutorial/ringdown.i
         block=Materials
         link=False
         language=python

shows the coefficients we will use. Note that $H_\mathrm{scale}$ is provided

A possible output of this tutorial problem using ParaView is provided below

!media media/ringdown_tut.png style=display:block;margin:auto;width:65%; caption=Left: Magnetic $\mathbf{m}$ (normalized) texture in glyphs and the color map showing the contrast of the demagnetizing potential $\Phi_\mathrm{H}$ at the conclusion of the ringdown. The vacuum block is shown as slightly transparent. Top Right: $\mathbf{m}$ components as a function of time during the ringdown. Bottom Right: Total energy as a function of time. id=mag_tutorial

!content pagination previous=tutorials/FE_phase_field_multi_domain_coupled.md next=tutorials/AFM_ringdown.md
