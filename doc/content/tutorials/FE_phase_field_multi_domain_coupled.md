!include tutorials/tutorial_header.md

# Tutorial 3: Ferroelectric thin film

This tutorial (and others) covers the basic usage of the ferroelectric (FE) phase field method implemented in FERRET. This specific example focuses on the thin film problem for $\mathrm{PbTiO}_3$ (PTO) or $\mathrm{BaTiO}_3$ (BTO) where periodicity is enforced along $x,y$ and $z$ corresponds to the film/substrate interface plane normal. We consider the calculation to be at room temperature.

As in Tutorial 2, we use the fully-coupled problem with electrostrictive free energy density. The total free energy density for this simulation can be written as,

\begin{equation}
  \begin{aligned}
    f = f_\mathrm{bulk} + f_\mathrm{elastic} + f_{\nabla P} + f_\mathrm{electostr} - \mathbf{P}\cdot\mathbf{E}
  \end{aligned}
\end{equation}

for the bulk, elastic, gradient, electrostrictive, and electrostatic free energies respectively. We choose a computational geometry as follows,

!listing tutorial/film.i
         block=Mesh
         link=False
         language=python

where $n_x, n_y,$ and $n_z$ are chosen accordingly such that the mesh spacing $\Delta = 1.0$ nm. In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through `Materials`, `Kernels`, or other MOOSE objects. For this problem, the length scale is introduced through the units in the `Materials` objects that connect to the `Kernels`. Note that this discretization is around the upper bound for these types of calculations since the thickness of the PTO wall is about $\Delta$. We use `SubdomainBoundingBoxGenerator` and `SideSetsBetweenSubdomainsGenerator` objects from MOOSE to split the rectilinear computational volume into two `block`s where the film sits on top of the substrate. For this example input file, we let the FE film have a thickness of 20 nm. This is a canonical problem in FE phase field problems (see [!cite](Li2001)) as it allows one to predict the domain topology as a function of temperature and film thickness.

!alert note title=Important!
We should note that the epitaxial strain coupling is not yet implemented in FERRET but this capability will be introduced in the future. As such, we just consider a free standing stress-free film.

The `Variables` block sets the `ICs` for $\langle \mathbf{P} \rangle = 0$ where the amplitude is set to be random fluctuations near 0

!listing tutorial/film.i
         block=Variables
         link=False
         language=python

Other variables $\Phi_\mathrm{E}, u_x, u_y$, and $u_z$ are also solved for. This tutorial problem evolves the time-dependent Landau-Ginzburg-Devonshire equation (TDLGD),

\begin{equation}
  \begin{aligned}
    \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta f}{\delta \mathbf{P}},
  \end{aligned}
\end{equation}

to find the ground state with $\Gamma_P = 1$ (arbitrary time scale). We also solve (at every time step) the conditions for electrostatic (Poisson equation) and mechanical equilibrium (stress divergence),

\begin{equation}
  \begin{aligned}
    \nabla \cdot \epsilon_b \nabla \Phi_\mathrm{E} &= \nabla \cdot \mathbf{P},\\
    \frac{\partial \sigma_{ij}}{\partial x_j} &= 0.
  \end{aligned}
\end{equation}

The variational derivatives of the total free energy density yield residual and jacobian contributions that are computed within the `Kernels` block,

!listing tutorial/film.i
         block=Kernels
         link=False
         language=python

We refer the reader to the following hyperlinks,

TDLGD:

- [`TimeDerivativeScaled`](source/kernels/TimeDerivativeScaled.md)
- [`BulkEnergyDerivativeEighth`](source/kernels/BulkEnergyDerivativeEighth.md)
- [`WallEnergyDerivative`](source/kernels/WallEnergyDerivative.md)
- [`Wall2EnergyDerivative`](source/kernels/Wall2EnergyDerivative.md)
- [`ElectrostrictiveCouplingPolarDerivative`](source/kernels/ElectrostrictiveCouplingPolarDerivative.md)
- [`PolarElectricPStrong`](source/kernels/PolarElectricPStrong.md)

Poisson equation:

- [`Electrostatics`](source/kernels/Electrostatics.md)
- [`PolarElectricEStrong`](source/kernels/PolarElectricEStrong.md)

Mechanical equilibrium:

- `TensorMechanics` (MOOSE `Action`) for $\partial \sigma_{ij} / \partial x_j = 0$
- [`ElectrostrictiveCouplingDispDerivative`](source/kernels/ElectrostrictiveCouplingDispDerivative.md)

for the different objects. Next, we define the boundary conditions on $\mathbf{P}$, $\mathbf{u}$ and $\Phi_\mathrm{E}$,

!listing tutorial/film.i
         block=BCs
         link=False
         language=python

where $x$ and $y$ are the lateral directions of the FE film. We utilize the [`GlobalStrain`](https://mooseframework.inl.gov/syntax/Modules/TensorMechanics/GlobalStrain/) system implemented in MOOSE to ensure periodicity of the strain tensor components along the periodic boundaries. This introduces a `ScalarKernel`,

\begin{equation}
  \begin{aligned}
    \int\limits_\Omega d^3 \mathbf{r} \sigma_{ij} = 0
  \end{aligned}
\end{equation}

with $\Omega$ the computational volume. We find a set of global displacement vectors `disp_x, disp_y, disp_z` (see `AuxKernels`) such that the above condition is satisfied (see [!cite](Biswas2020) for more description of the method). Note that the substrate has a boundary condition such that the elastic displacement vector $\mathbf{u}$ goes to zero deep in the substrate which is handled by the `DirichletBC`. Some possible outputs of this problem with PTO `Materials` coefficients at room temperature,

!listing tutorial/film.i
         block=Materials
         link=False
         language=python

are visualized below using ParaView Filters for colormaps of $|\mathbf{P}|$ and white arrow Glyphs for the $\mathbf{P}$ directors.

!media media/tut_film.png style=display:block;margin:auto;width:50%; caption=$|\mathbf{P}|$ across the PTO film sample in 3D. The Glyph filter provides the arrows corresponding to the polar director.   id=fig-ferret_tut3

If we use the BTO `Materials` coefficients at room temperature, we have the below output.

!media media/tut_film2.png style=display:block;margin:auto;width:50%; caption=$|\mathbf{P}|$ across the BTO film sample in 3D. The Glyph filter provides the arrows corresponding to the polar director.   id=fig-ferret_tut4

In each of the calculations, the boundary condition on $\Phi_\mathrm{E}$ is open-circuit (free) at the surface of the thin film. This leads to a strong electrostatic action to form in-plane domains. Both input files are provided in the tutorial folder in the FERRET root directory.

!content pagination previous=tutorials/ferroelectric_domain_wall.md next=tutorials/magnetic_ringdown.md
