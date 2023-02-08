!include tutorials/tutorial_header.md

# Tutorial 2: Ferroelectric domain wall

This tutorial gives an example on how to compute the $180^\circ$ domain wall (DW) profile of $\mathrm{BaTiO}_3$ (BTO) at room temperature ($T = 298$ K) in FERRET. In this problem, the polarization $\mathbf{P}$ is coupled to the elastic strain field $\varepsilon_{kl} = \frac{1}{2} \left( \frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k} \right)$ via electrostrictive coupling appropriate for BTO.

The expected DW plane configuration is (100) or equivalent directions so we choose a computational geometry as follows,



where $n_x, n_y,$ and $n_z$ are chosen accordingly such that the mesh spacing $\Delta$ is approximately 0.25 nm.


We utilize the GlobalStrain system in MOOSE to ensure periodicity of the strain tensor components along the long direction of the box.


In principle, this type of calculation can be generalized to any ferroic material (i.e. ferromagnets) to study the DW textures of magnetization in the presence of additional couplings (for example magnetoelasticity).



run time 1864 secs on 8 processors

TEST: 
