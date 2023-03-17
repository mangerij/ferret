!config navigation breadcrumbs=False

# H2020-MSCA-IF-2019, Example: ground states

!alert construction title=Documentation in-progress
This section will have information that will be available when the review process is complete for the main publication encompassing this work. Please contact the developers if you need assistance with this aspect of the module.

This page details how to obtain homogeneous solutions of BFO for the structural and spin order. The first calculation (`BFO_homogeneous_PA.i`) determines the polarization $\mathbf{P}$ and antiphase octahedral cage tilt vectors $\mathbf{A}$ along with the coupled elastic strain tensor components $\varepsilon_{ij}$ computed from the elastic displacement variable $\mathbf{u}$.

The second (`BFO_homogeneous_LM.i`) uses information from the first to determine the antiferromagnetic order corresponding to the Neel vector $\mathbf{L}$ and total magnetization $\mathbf{m}$. We can separate these calculations because the magnetic system does not appreciably influence the relaxation dynamics of the polarization since the work required to move ions (i.e. those relative displacements corresponding to $\mathbf{P}$ and $\mathbf{A}$) is much larger than those to reorient the Fe spins.

# Homogeneous structural order

We first consider a simple $2\times 2\times 2$ nm on a grid of $5 \times 5 \times 5$ finite elements with the `Mesh` block,

!listing tutorial/BFO_homogeneous_PA.i
         block=Mesh
         link=False
         language=python

The `Variables` block defines the relevant variables $\{\mathbf{P}, \mathbf{A}, \mathbf{u}\}$ of the simulation,

!listing tutorial/BFO_homogeneous_PA.i
         block=Variables
         link=False
         language=python

We also include global strain quantities defined by the relationship,

\begin{equation}
  \begin{aligned}
    \varepsilon_{ij} = \varepsilon_{ij}^\mathrm{global} + \varepsilon_{ij}^\mathrm{local}.
  \end{aligned}
\end{equation}

In this problem, the local strain is zero (no inhomogeneities) which means that the spontaneous strain tensor arising from the electro- and rotostrictive coupling is equal to the global strain tensor. We utilize the `GlobalStrain` system implemented in MOOSE to ensure periodicity of the total strain tensor components (global+local) along the periodic boundaries of the box ($i = x,y,z, j = x,y,z$). This introduces a `ScalarKernel`,

\begin{equation}
  \begin{aligned}
    \int\limits_\Omega d^3 \mathbf{r} \sigma_{ij} = 0
  \end{aligned}
\end{equation}

with $\Omega$ the computational volume. We find a set of global displacement vectors `disp_x, disp_y, disp_z` (see `AuxKernels`) such that the above condition is satisfied (see [!cite](Biswas2020) for an extended description of the method). In the input file, we see a number of objects that allow for this additional system,

!listing tutorial/BFO_homogeneous_PA.i
         block=ScalarKernels
         link=False
         language=python

and

!listing tutorial/BFO_homogeneous_PA.i
         block=UserObjects
         link=False
         language=python

To find the lowest energy solution (the ground state), we evolve the the time dependent Landau-Ginzburg equations,

\begin{equation}\label{eqn:TDLG}
  \begin{aligned}
    \frac{\partial \mathbf{P} }{\partial t} &= -\Gamma_P \frac{\delta f}{\delta \mathbf{P}}, \\
    \frac{\partial \mathbf{A} }{\partial t} &= -\Gamma_A \frac{\delta f}{\delta \mathbf{A}},
  \end{aligned}
\end{equation}

along with the equation for mechanical equilibrium for the elastic strain field $\varepsilon_{kl}$,

\begin{equation}
  \begin{aligned}
    \frac{\partial \sigma_{ij}}{\partial x_j} = \frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + Q_{klmn} P_m P_n + R_{klmn} A_m A_n \right)\right] = 0,
  \end{aligned}
\end{equation}

where $Q_{ijkl}$ and $R_{ijkl}$ are the electrostrictive and rotostrictive coefficient tensors with $\mathbf{P}$ and $\mathbf{A}$ the order parameters associated with the spontaneous polarization and antiphase tilting of the oxygen octahedral cages. The mechanical equilibrium equation is assumed to be satisfied for every time step which is a reasonable approximation for the elastic strains that arise during domain evolution. The total free energy density, $f$ corresponds to the fourth-order coupled potential parameterized by the DFT work detailed [!cite](Fedorova2022). We set $\Gamma_P = \Gamma_A$ to unity for these calculations since we are only interested in the final state. The `Kernels` block,

!listing tutorial/BFO_homogeneous_PA.i
         block=Kernels
         link=False
         language=python

sets up the necessary residuals and jacobians for these equations. The initial conditions are set such that $\mathbf{P}$ and $\mathbf{A}$ are close to the expected values.


 and the boundary conditions (periodic) are both chosen to produce a homogeneous solution with $\mathbf{P}$ and $\mathbf{A}$ parallel ($\mathbf{P}\uparrow\uparrow\mathbf{A}$). The magnitudes of the structural order parameters are $P_s = |\mathbf{P}| = 0.945 \,\,\mathrm{C}/\mathrm{m}^2$  and $A_s = |\mathbf{A}| = 13.398^\circ$. The spontaneous $homogeneous$ normal and shear strain values are listed as rows with $\varepsilon_n = 1.308\times10^{-2}$ and $\varepsilon_s = 2.95 \times 10^{-3}$ respectively. The eight-fold domain variant symmetry is found in the pseudocubic reference frame corresponding to the below table.

!equation
\begin{array}{c|cccccccc}
\hline
\hline
 & & & & & & & & \\
 (\mathbf{P}\uparrow\uparrow\mathbf{A}): & [111] & [\bar{1}11] & [1\bar{1}\bar{1}] & [\bar{1}\bar{1}\bar{1}] & [1\bar{1}1] & [11\bar{1}] & [\bar{1}\bar{1}1] & [\bar{1}1\bar{1}] \\
 \varepsilon_{xx}: & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n\\
 \varepsilon_{yy}: & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n\\
 \varepsilon_{zz}: & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n & \varepsilon_n\\
 \varepsilon_{xy}: & \varepsilon_s & -\varepsilon_s & -\varepsilon_s & \varepsilon_s & -\varepsilon_s & \varepsilon_s & \varepsilon_s & -\varepsilon_s \\
 \varepsilon_{yz}: & \varepsilon_s & \varepsilon_s & \varepsilon_s & \varepsilon_s &-\varepsilon_s & -\varepsilon_s & -\varepsilon_s & -\varepsilon_s \\
 \varepsilon_{xz}: & \varepsilon_s & -\varepsilon_s & -\varepsilon_s & \varepsilon_s & \varepsilon_s & -\varepsilon_s & -\varepsilon_s & \varepsilon_s \\
\end{array}

Note that $(\mathbf{P} \uparrow \downarrow \mathbf{A})$ is also a possible minimized energy solution of the thermodynamic potential $f$. Due to the symmetry of the electrostrictive and rotostrictive coupling terms, the table is left invariant under full reversal of $\mathbf{A}$. The free energy density of the eight-fold domain possibilities is -15.5653 $\mathrm{eV}\cdot\mathrm{nm}^{-3}$.


# Homogeneous spin order

For this example, please consult `BFO_P0A0_mRD.i` in the tutorials subdirectory. We consider the output of the previous simulation $\{\mathbf{P},\mathbf{A}\}$ as an initial condition to find the homogeneous spin ground state. To do this, we use the following `Mesh` and `Variables` blocks, 

Our variables to solve for are the set of sublattice magnetizations $\{\mathbf{m}_1, \mathbf{m}_2\}$ with postprocessed outputs $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$. Note that the conventional factor of 2 is missing from the output of these calculations.



This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_results2.md next=MSCA_EU_Horizon2020_Results/horizon2020_ex2.md
