!include tutorials/tutorial_header.md

# Tutorial 5: Piezoelectric actuation

This tutorial covers the basic usage of the piezoelectric coupling implemented in FERRET.

In the linear limit, the constitutive governing equation for piezoelectricity is

\begin{equation}
  \begin{aligned}
    \frac{\partial}{\partial x_j} \left(C_{ijkl} \varepsilon_{kl} + d_{klm} E_m\right) = 0,
  \end{aligned}
\end{equation}

where $C_{ijkl}, \varepsilon_{kl},$ and $E_m$ are the components of the elastic stiffness tensor, elastic strain tensor, and electric field. The direct piezoeletric coefficient $d_{klm}$ is of rank three. This is equivalently the equation for mechanical equilibrium with coupling to the electric field,  $E_m = - \nabla \Phi_\mathrm{E}$. Note that the strain tensor is defined in the usual way, $\varepsilon_{kl} = \frac{1}{2} \left( \frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k} \right).$

Additionally, the Poisson equation also includes contributions from the (converse piezo) strain-charge,

\begin{equation}
  \begin{aligned}
    \frac{\partial }{\partial x_r } \epsilon_b \frac{\partial \Phi_\mathrm{E}}{\partial x_r} + \frac{\partial d_{jkl} \sigma_{kl}}{\partial x_j} = 0,
  \end{aligned}
\end{equation}

 with $\sigma_{kl}$ being the components of the elastic stress tensor. With sufficient choice of materials parameters, Equations (1) and (2) can be solved self-consistently under arbitrary mechanical loads or applied electric fields to yield the static configuration of the electrostatic potential and elastic displacements $\mathbf{u}.$

These equations can be cast dynamically, to simulate piezeoelectric actuation in real time. This is done by setting the RHS of the mechanical equilibrium equation to be equal to $\partial u_i / \partial t$

In this problem, we consider a computational geometry

!listing tutorial/piezoelectric.i
         block=Mesh
         link=False
         language=python

 The finite element mesh discretization schema is chosen to be quadrilateral elements HEX8. In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through Materials, Kernels, or other MOOSE objects. We seed the materials coefficients $C_{ijkl}$ and $d_{ijk}$ through the following block.

!listing tutorial/piezoelectric.i
         block=Materials
         link=False
         language=python

which are given in units of GPa (for $C_{ijkl}$). Also listed is the permitivitty of the medium denoted as $\epsilon_b$ and some Materials objects to compute the linear elastic strain and stress. The Kernels suitable for solving the above set of equations are

!listing tutorial/piezoelectric.i
         block=Kernels
         link=False
         language=python

which are for the constitutive governing mechanical equations of the piezoelectric (the TensorMechanics Action which sets up $\partial \sigma_{ij} / \partial x_j$ and ConversePiezoelectricStrain which handles the coupling). Also in the Kernels block, we have the Poisson equation, Electrostatics (\nabla^2 \Phi_\mathrm{E}) and PiezoelectricStrainCharge which handles the coupling to the bound charge arising from the strain field.

Finally we have some optional AuxVariables that are computed in

!listing tutorial/piezoelectric.i
         block=AuxsKernels
         link=False
         language=python

which are stored to be viewed in the output. Some possible outputs of this tutorial problem could look like the figure below.

!media media/piezo_tutorial.png style=display:block;margin:auto;width:50%; caption=Top: Warped (x100) Filter showing the modulation of the displacement vectors under the applied electric field. Bottom: $\sigma_{xx}$ and $\sigma_{zz} components as a function of time during the actuation. id=piezo_tutorial
