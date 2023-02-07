!include tutorials/tutorial_header.md

# Tutorial 5: Piezoelectric actuation

This tutorial covers the basic usage of the piezoelectric coupling implemented in FERRET.

In the linear limit, the constitutive governing equation for piezoelectricity is

\begin{equation}
  \begin{aligned}
    \frac{\partial}{\partial x_j} \left(C_{ijkl} \varepsilon_{kl} + d_{klm} E_m\right) = 0,
  \end{aligned}
\end{equation}

where $C_{ijkl}, \varepsilon_{kl},$ and $E_m$ are the components of the elastic stiffness tensor, elastic strain tensor, and electric field. The direct piezoeletric coefficient $d_{klm}$ is of rank three. Note that the strain tensor is defined in the usual way, $\varepsilon_{kl} = \frac{1}{2} \left( \frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k} \right).$

Additionally, the Poisson equation also includes contributions from the (converse piezo) strain-charge,

\begin{equation}
  \begin{aligned}
    \frac{\partial^2 \Phi_\mathrm{E} }{\partial x_r \partial x_r} + \frac{\partial d_{jkl} \sigma_{kl}}{\partial x_j} = 0 ,
  \end{aligned}
\end{equation}

 with $\sigma_{kl}$ being the components of the elastic stress tensor. With sufficient choice of materials parameters, Equations (1) and (2) can be solved self-consistently under arbitrary mechanical loads or applied electric fields to yield the static configuration of the electrostatic potential and elastic displacements $\mathbf{u}.$

Equations (1) and (2) can be cast dynamically, to simulate piezeoelectric actuation in real time.

In this problem, we consider a computational geometry

!listing tutorial/piezoelectric.i
         block=Mesh
         link=False
         language=python

In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through Materials, Kernels, or other MOOSE objects. We seed the materials coefficients $C_{ijkl}$ and $r_{ijk}$ through the following block.

!listing tutorial/piezoelectric.i
         block=Materials
         link=False
         language=python
