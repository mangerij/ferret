!include tutorials/tutorial_header.md

# Tutorial 1: Multidomain state of a FE block

This tutorial (and others) covers the basic usage of the ferroelectric phase field method implemented in FERRET.

This simple problem considers solving the time-dependent LGD equation with just bulk and gradient contributions to the free energy in the presence of a depolarization phenomena (via solving the Poisson equation iteratively).

Consider a computational domain with a geometry $(30\times 30\times 6)$ which we define with the Mesh block.

!listing tutorial/multidomain.i
         block=Mesh
         link=False
         language=python

In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through `Materials`, `Kernels`, or other MOOSE objects. For this problem, the length scale is introduced through the units in the `Materials` objects that connect to the `Kernels`. Ferret uses a special base units system for a number of problems. This reduces the load quite extensively on the PETSc solvers to iterate the problem.

In the total free energy density, we include the bulk free energy density,

\begin{equation}\label{bulk}
\begin{aligned}
f_\mathrm{bulk} &= \alpha_1 \left(P_x^2 + P_y^2 + P_z^2 \right) + \alpha_{11} \left(P_x^4 + P_y^4 + P_z^4 \right)  + \alpha_{12} \left(P_x^2 P_y^2 + P_y^2 P_z^2 + P_x^2 P_z^2 \right) + \alpha_{111} \left(P_x^6 + P_y^6 + P_z^6 \right) \\
&+ \alpha_{112} \left[P_x^4 \left(P_y^2 + P_z^2 \right) + P_y^4 \left(P_x^2 + P_z^2 \right) + P_z^4 \left(P_x^2 + P_y^2 \right) \right] + \alpha_{123} \left(P_x^2 P_y^2 P_z^2 \right). \\
\end{aligned}
\end{equation}

to the sixth order and the wall free energy density, $f_\mathrm{wall}$,

\begin{equation}\label{wall}
\begin{aligned}
f_\mathrm{wall} &= \frac{1}{2} G_{11} \left[\left(\frac{\partial P_x}{\partial x} \right)^2  +  \left(\frac{\partial P_y}{\partial y}\right)^2  + \left( \frac{\partial P_z}{\partial z} \right)^2\right] + G_{12} \left[\frac{\partial P_x}{\partial x}\frac{\partial P_y}{\partial y} +  \frac{\partial P_y}{\partial y} \frac{\partial P_z}{\partial z} +  \frac{\partial P_x}{\partial x} \frac{\partial P_z}{\partial z} \right] \\
&+ \frac{1}{2} G_{44} \left[\left(\frac{\partial P_x}{\partial y}  +   \frac{\partial P_y}{\partial x} \right)^2 +\left(\frac{\partial P_y}{\partial z}  +  \frac{\partial P_z}{\partial y} \right)^2 + \left(\frac{\partial P_x}{\partial z}  +  \frac{\partial P_z}{\partial x} \right)^2  \right] \\
&+ \frac{1}{2} G'_{44} \left[\left(\frac{\partial P_x}{\partial y}  -  \frac{\partial P_y}{\partial x} \right)^2  +  \left(\frac{\partial P_y}{\partial z}  -  \frac{\partial P_z}{\partial y} \right)^2 +  \left(\frac{\partial P_x}{\partial z}  -  \frac{\partial P_z}{\partial x} \right)^2 \right] \\
\end{aligned}
\end{equation}

where $\alpha_{ijk}$ and $G_{ij}$ are suitable coefficients for lead-titanate. To solve the problem, we need the time-dependent Landau-Ginzburg equation,

\begin{equation}\label{tdlgd}
\begin{aligned}
  \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta F}{\delta \mathbf{P}},
\end{aligned}
\end{equation}

which has variational derivatives of the total free energy,

\begin{equation}\label{varderiv}
\begin{aligned}
  \frac{\delta F}{\delta \mathbf{P}} = \frac{\partial F}{\partial \mathbf{P}} - \frac{\partial}{\partial x_k} \left(\frac{\partial F}{\partial \left(\frac{\partial \mathbf{P}}{\partial x_k}\right)} \right)
\end{aligned}
\end{equation}

The computations of the two relevant energy terms are introduced through a number of `Kernels` in the input files as below.

!listing tutorial/multidomain.i
         block=Kernels
         link=False
         language=python

We need one Kernel for each variable ($P_x, P_y, P_z$, and $\Phi_\mathrm{E}$). The computation of Eq. \ref{varderiv} is split into three Kernel contributions for each of the three components of $\mathbf{P}$. The FERRET Syntax page details the algebra needed to turn these terms ([`BulkEnergyDerivativeEighth`](source/kernels/BulkEnergyDerivativeEighth.md), [`WallEnergyDerivative`](source/kernels/WallEnergyDerivative.md), and [`Wall2EnergyDerivative`](source/kernels/Wall2EnergyDerivative.md)) into weak forms.

The reader will notice that we also have included terms related to the free energy due to an applied field $-\mathbf{P}\cdot\mathbf{E}$ in [`PolarElectricPStrong`](source/kernels/PolarElectricPStrong.md) and also terms related to the coupling to the electric field. The Poisson equation for our system reads,

\begin{equation}\label{poisson}
\begin{aligned}
 \nabla \cdot \epsilon_b \nabla \Phi_\mathrm{E} = \nabla \cdot \mathbf{P}
\end{aligned}
\end{equation}

which is satisfied at every step of the time evolution. The left hand side is computed by Electrostatics and the right hand side computed by [`PolarElectricPStrong`](source/kernels/PolarElectricPStrong.md). Note that $\epsilon_b$ corresponds to a background dielectric strength which is assigned to high frequency (core) electrons. Next, we have the time derivatives of $P_x, P_y$, and $P_z$ with [`TimeDerivativeScaled`](source/kernels/TimeDerivativeScaled.md). We implement some reduction of input file length by using the `GlobalParams` block.

!listing tutorial/multidomain.i
         block=GlobalParams
         link=False
         language=python

which adds these lines automatically in the relevant `Kernel` and `Materials` objects. To seed the relevant materials coefficients to the problem, we use the Materials system in MOOSE,

!listing tutorial/multidomain.i
         block=Materials
         link=False
         language=python

In this problem, the boundary conditions are left empty. This is equivalently known as a natural boundary condition in MOOSe. As such, large depolarizing fields are generated from the surface charges ($\sigma_b = \mathbf{P}\cdot \mathbf{n}$) where $\mathbf{n}$ is some surface normal. The action of the depolarization drives the polar topology into a classic flux-closure domain pattern as shown in the below figure.

!media media/ferret-tut1.png style=display:block;margin:auto;width:50%; caption=Classic flux-closure domain pattern solution to the first tutorial.  id=fig-ferret_tut1

A few other FERRET and MOOSE objects exist in the input file. We use the `Postprocessors` system to track various aspects of the simulation.

!listing tutorial/multidomain.i
         block=Postprocessors
         link=False
         language=python

The bulk and gradient energy contributions are listed (computed volume integral over the simulation box). Their sum is computed with the LinearCombinationPostprocessor to be equal to $F_\mathrm{tot}$. An astute user should ask why the electric field contribution $-\mathbf{P}\cdot\mathbf{E}$ is also not included here. We leave that as an exercise for the user to verify that the total energy with this contribution still is decreasing with every time step as gauranteed by the mathematical form of Eq. \refLGD.

Then, the relative change of this energy between adjacent time steps is calculated with PercentChangePostprocessor as,

\begin{equation}\label{rel}
\begin{aligned}
  r = \frac{F_\mathrm{tot}(t_k) - F_\mathrm{tot}(t_{k-1})}{dt}
\end{aligned}
\end{equation}

The Terminator located in the UserObjects block

!listing tutorial/multidomain.i
         block=UserObjects
         link=False
         language=python

which tracks when the postprocessed value $r$ is less than $10^{-4}$. This is a useful tool to end the simulation when the suspected ground state has been reached. In general, this value $10^{-8} < r < 10^{-4}$ seems to be sufficient for these types of problems.

Finally, we share our PETSc and `Executioner` options that seem to be the most efficient for polar domain prediction problems.

!listing tutorial/multidomain.i
         block=Preconditioning
         link=False
         language=python

and

!listing tutorial/multidomain.i
         block=Executioner
         link=False
         language=python

These flags can be optimized to lead to less wall clock times but this is a subject of future work.

!content pagination next=tutorials/ferroelectric_domain_wall.md
