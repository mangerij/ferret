!include tutorials/tutorial_header.md

# Tutorial 2: Ferroelectric domain wall

This tutorial gives an example on how to compute the $180^\circ$ domain wall (DW) profile of $\mathrm{BaTiO}_3$ (BTO) at room temperature ($T = 298$ K) in FERRET. In this problem (as opposed to Tutorial 1), the polarization $\mathbf{P}$ is coupled to the elastic strain field $\varepsilon_{kl} = \frac{1}{2} \left( \frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k} \right)$ via electrostrictive coupling appropriate for BTO. The variable $u_k$ is the $k^\mathrm{th}$ component of the elastic displacement vector $\mathbf{u}$. The total free energy density is

\begin{equation}
  \begin{aligned}
    f = f_\mathrm{bulk} + f_\mathrm{elastic} + f_{\nabla P} + f_\mathrm{electostr} - \mathbf{P}\cdot\mathbf{E}
  \end{aligned}
\end{equation}

for the bulk, elastic, gradient, electrostrictive, and electrostatic free energies respectively. The expected DW plane configuration is (100) or equivalent directions so we choose a computational geometry as follows,

!listing tutorial/ferroelectric_domain_wall.i
         block=Mesh
         link=False
         language=python

where $n_x, n_y,$ and $n_z$ are chosen accordingly such that the mesh spacing $\Delta = 0.5$ nm. In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through `Materials`, `Kernels`, or other MOOSE objects. For this problem, the length scale is introduced through the units in the `Materials` objects that connect to the `Kernels`. Ferret uses a special base units system (aC, kg, nm, sec) for a number of problems related to ground state prediction. Note that for BTO, $P_s = 0.26\,\mathrm{C}/\mathrm{m}^2 = 0.26 \,\mathrm{aC}/\mathrm{nm}^2$. This unit system choice reduces the load quite extensively on the PETSc solvers to iterate the problem. To simulate the DW texture, we choose a $\cos(x)$ profile for $P_z$ with the `ICs` block inside the `Variables` block as,

!listing tutorial/ferroelectric_domain_wall.i
         block=Variables
         link=False
         language=python

where the `FunctionIC` defines a function called `DW_func` with,

\begin{equation}
  \begin{aligned}
    P_z = 0.2 \,\cos{\left(\frac{2\pi}{L_x} x\right)}.
  \end{aligned}
\end{equation}

Other variables $P_x, P_y, \Phi_\mathrm{E}, u_x, u_y$, and $u_z$ are also solved for. This tutorial problem evolves the time-dependent Landau-Ginzburg-Devonshire equation (TDLGD),

\begin{equation}
  \begin{aligned}
    \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta f}{\delta \mathbf{P}},
  \end{aligned}
\end{equation}

to find the ground state with $\Gamma_P = 1$ (arbitrary time scale). We also solve (at every time step) the conditions for electrostatic (Poisson equation) and mechanical equilibrium (stress divergence),

\begin{equation}
  \begin{aligned}
    \nabla \cdot \epsilon_b \nabla \Phi_\mathrm{E} &= \nabla \cdot \mathbf{P},\\
    \frac{\partial \sigma_{ij}}{\partial x_j} &= 0,
  \end{aligned}
\end{equation}

with $\epsilon_b$ a background dielectric permittivity. The variational derivatives of the total free energy density yield residual and jacobian contributions that are computed within the `Kernels` block,

!listing tutorial/ferroelectric_domain_wall.i
         block=Kernels
         link=False
         language=python

The weak-form algebra required for setting up these objects are provided in the FERRET syntax list. An interested user may click the hyperlinks here:

TDLGD:

- [`TimeDerivativeScaled`](source/kernels/TimeDerivativeScaled.md)
- [`BulkEnergyDerivativeEighth`](source/kernels/BulkEnergyDerivativeEighth.md)
- [`WallEnergyDerivative`](source/kernels/WallEnergyDerivative.md)
- [`CubicParentElasticPDerivative`](source/kernels/CubicParentElasticPDerivative.md)
- [`PolarElectricPStrong`](source/kernels/PolarElectricPStrong.md)

Poisson equation:

- [`Electrostatics`](source/kernels/Electrostatics.md)
- [`PolarElectricEStrong`](source/kernels/PolarElectricEStrong.md)

Mechanical equilibrium:

- `SolidMechanics` (MOOSE `Action`) for $\partial \sigma_{ij} / \partial x_j = 0$

for the different objects. The electrostrictive coupling to the mechanics enters through the `Materials` block as an eigenstrain rather than through a separate coupling kernel, so no `ElectrostrictiveCoupling*` kernel appears in this input file. Note that the `Materials` block via

!listing tutorial/ferroelectric_domain_wall.i
         block=Materials
         link=False
         language=python

We utilize the `GlobalStrain` system implemented in MOOSE to ensure periodicity of the strain tensor components along the long direction of the box ($i = x, j = x$). This introduces a `ScalarKernel`,

\begin{equation}
  \begin{aligned}
    \int\limits_\Omega d^3 \mathbf{r} \sigma_{ij} = 0
  \end{aligned}
\end{equation}

with $\Omega$ the computational volume. We find a set of global displacement vectors `disp_x, disp_y, disp_z` (see `AuxKernels`) such that the above condition is satisfied (see [!cite](Biswas2020) for an extended description of the method). In the input file, these objects are generated by the `GlobalStrain` action,

!listing tutorial/ferroelectric_domain_wall.i
         block=Physics/SolidMechanics/GlobalStrain
         link=False
         language=python

The action creates the `ScalarKernel`, its associated user object and the auxiliary displacement fields automatically, so no explicit `ScalarKernels` or `UserObjects` entries are required for the global strain system. It works together with the `BCs` block,

!listing tutorial/ferroelectric_domain_wall.i
         block=BCs
         link=False
         language=python

which ensures the appropriate periodicity along the long direction of the box. We find that setting periodicity along the $y$ and $z$ directions does not influence the system variables and is a redundant boundary condition. This problem also has a number of `AuxVariables` to store the elastic strain and global displacement fields. Finally, it should be noted that in the `UserObjects` block, we also include a `Terminator` object which kills the problem when the relative change of the total energy between adjacent time steps is less than $5\times 10^{-7}$. The wall clock time of this problem is 316.8 secs on 6 processors. After the problem is solved, a typical output can be viewed in ParaView as below.

!media media/DW_prof.png style=display:block;margin:auto;width:50%; caption=Top: $P_z$ across the DW region. Bottom: Variation of $\varepsilon_{xx}$ and $\varepsilon_{yy}$ along the same arclength ($x$) in nanometers.   id=fig-ferret_tut2

showing the thickness of the DW region along with the variations of the spontaneous strains $\varepsilon_{xx}$ and $\varepsilon_{yy}$. The resulting order parameters in the homogeneous region are $P_s = P_z = 0.26523$ $\mathrm{C}/\mathrm{m}^2$, and normal strains $\varepsilon_{zz} = 7.6485\times 10^{3}$ and $\varepsilon_{xx} = \varepsilon_{yy} = -3.12893 \times 10^{-3}$ which is in good agreement with the results of [!cite](Hlinka2006). The thickness, which can be calculated by fitting $P_z$ to a $tanh(x)$ profile, agrees well with the calculations from [!cite](Marton2010) which highlights a number of a different BTO DWs as a function of temperature.

In principle, this type of calculation can be generalized to any ferroic material (i.e. ferromagnets or multiferroics) to study the DW textures of order parameters in the presence of additional couplings (for example magnetoelasticity or the flexoelectric coupling to gradients in the strain field).

!content pagination previous=tutorials/FE_phase_field_multi_domain.md next=tutorials/FE_phase_field_multi_domain_coupled.md
