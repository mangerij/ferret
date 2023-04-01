!config navigation breadcrumbs=False

# H2020-MSCA-IF-2020, Example: ground states

This page details how to obtain homogeneous solutions of BFO for the structural and spin order. The first calculation (`BFO_homogeneous_PA.i`) determines the polarization $\mathbf{P}$ and antiphase octahedral cage tilt vectors $\mathbf{A}$ along with the coupled elastic strain tensor components $\varepsilon_{ij}$ computed from the elastic displacement variable $\mathbf{u}$.

The second (`BFO_P0A0_mRD.i`) uses information from the first to determine the antiferromagnetic order corresponding to the Neel vector $\mathbf{L}$ and total magnetization $\mathbf{m}$. We can separate these calculations because the magnetic system does not appreciably influence the relaxation dynamics of the polarization since the work required to move ions (i.e. those relative displacements corresponding to $\mathbf{P}$ and $\mathbf{A}$) is much larger than those to reorient the Fe spins.

# Homogeneous structural order

We first consider a simple $1\times 1\times 1$ nm on a grid of $3 \times 3 \times 3$ finite elements with the `Mesh` block,

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

with $\Omega$ the computational volume. We find a set of global displacement vectors `disp_x, disp_y, disp_z` (see `AuxKernels`) such that the above condition is satisfied (see [!cite](Biswas2020) for an extended description of the method). In the input file, we see a number of MOOSE objects that allow for this additional system,

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

sets up the necessary residuals and jacobians for these equations. The initial conditions are set such that $\mathbf{P}$ and $\mathbf{A}$ are close to the expected values. The boundary conditions (periodic) are both chosen to produce a homogeneous solution with $\mathbf{P}$ and $\mathbf{A}$ parallel ($\mathbf{P}\uparrow\uparrow\mathbf{A}$). We set this condition with the `BCs` block,

!listing tutorial/BFO_homogeneous_PA.i
         block=BCs
         link=False
         language=python

Note that the `DirichletBC` condition $\mathbf{u} = 0$ is set so there are no translations as $\mathbf{P}$ and $\mathbf{A}$ evolve to their minimum energy values. The resulting magnitudes of the structural order parameters are $P_s = |\mathbf{P}| = 0.945 \,\,\mathrm{C}/\mathrm{m}^2$  and $A_s = |\mathbf{A}| = 13.398^\circ$. The spontaneous $homogeneous$ normal and shear strain values are listed as rows with $\varepsilon_n = 1.308\times10^{-2}$ and $\varepsilon_s = 2.95 \times 10^{-3}$ respectively. The eight-fold domain variant symmetry is found in the pseudocubic reference frame corresponding to the below table.

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

Note that $(\mathbf{P} \uparrow \downarrow \mathbf{A})$ are also possible minimized energy solutions of the thermodynamic potential $f$. Due to the symmetry of the electrostrictive and rotostrictive coupling terms, the table is left invariant under full reversal of $\mathbf{A}$. The free energy density of the eight-fold domain possibilities is -15.5653 $\mathrm{eV}\cdot\mathrm{nm}^{-3}$.

This simulation runs in 18.23 seconds on 4 processors using the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) distribution of MOOSE.

# Homogeneous spin order

For this example, please consult `BFO_P0A0_mRD.i` in the tutorials subdirectory. We consider the output of the previous simulation $\{\mathbf{P},\mathbf{A}\}$ as an initial condition to find the homogeneous spin ground state. To do this, we use the following `Mesh` block,

!listing tutorial/BFO_P0A0_mRD.i
         block=Mesh
         link=False
         language=python

and `AuxVariables` blocks,

!listing tutorial/BFO_P0A0_mRD.i
         block=AuxVariables
         link=False
         language=python

where we use (i.e. for vector component $P_z$) `initial_from_file_var = polar_z` and `initial_from_file_timestep = 'LATEST'` as code flags to load in the values of `polar_z` and others at the latest timestep from the Exodus output `BFO_P0A0.e`. In this sense, we have changed $\{\mathbf{P},\mathbf{A}\}$ to be an AuxVariable class instead of a Variable since they will be fixed for the duration of the evolution of the spins. Our variables to solve for are the set of sublattice magnetizations $\{\mathbf{m}_1, \mathbf{m}_2\}$ with postprocessed outputs of the Néel vector $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and net (weak) magnetization $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$. Note that the conventional factor of 2 is missing from the output of these calculations.

This problem considers the evolution of the normalized Landau-Lifshitz-Gilbert (Bloch) equation for the two sublattices,

\begin{equation}\label{eqn:LLG_LLB}
  \begin{aligned}
    \frac{d \mathbf{m}_\eta}{dt} = -\frac{\gamma}{1+\alpha^2} \left(\mathbf{m}_\eta \times \mathbf{H}_\eta\right) - \frac{\gamma\alpha}{1+\alpha^2}\mathbf{m}_\eta\times\left(\mathbf{m}_\eta\times\mathbf{H}_\eta\right) + \frac{\gamma\tilde\alpha_\parallel}{(1+\alpha^2)} m_{\eta}^2 \left[ m_{\eta}^2 - 1 \right]\mathbf{m}_\eta,
  \end{aligned}
\end{equation}

where $\gamma$ corresponds to the electron gyromagnetic factor equal to $2.2101\times10^5$ rad. m A${}^{-1}$ s ${}^{-1}$, $\alpha$ the Gilbert damping, and $\tilde\alpha_\parallel$ the longitudinal damping constant from the LLB approximation. The effective fields are defined as $\mathbf{H}_\eta = - \mu_0^{-1} M_s^{-1} \delta f_\mathrm{mag} / \delta \mathbf{m}_\eta$ with $\mu_0$ the permeability of vacuum. The saturation magnetization density of the BFO sublattices is $M_s = 4.0$ $\mu$B/Fe [!cite](Dixit2015). We choose an initial condition close to the expected ground state values. We set $\alpha = 0.1$ and $\tilde\alpha_\parallel = 10^3$ and ringdown the magnetic system to an energy minimum.

The spin free energy, $f_\mathrm{mag}$, is correspondingly,

\begin{equation}\label{eqn:fmag}
  \begin{aligned}
    f_\mathrm{mag} = f_\mathrm{exch} + f_\mathrm{DMI} + f_\mathrm{easy} + f_\mathrm{anis}
  \end{aligned}
\end{equation}

where

\begin{equation}\label{eqn:fmag_split}
  \begin{aligned}
    f_\mathrm{exch} &= D_e \left(\mathbf{L}^2 - \mathbf{m}^2\right)\\
    f_\mathrm{DMI} &= 4 D_0 \mathbf{A} \cdot \left( \mathbf{L} \times \mathbf{m}\right) \\
    f_\mathrm{easy} &= \sum\limits_{\eta = 1}^2 K_1 \left(\textbf{m}_\eta \cdot \hat{\mathbf{P}}\right)^2 \\
    f_\mathrm{anis} &= \sum\limits_{\eta = 1}^2 \left(K_1^c + a|\mathbf{A}|^2\right)\left(m_{\eta,x}^2 m_{\eta,y}^2 m_{\eta,z}^2\right).
  \end{aligned}
\end{equation}

Each of these contributions to the effective fields $\mathbf{H}$ contain different physics. For example, the effective fields from $f_\mathrm{exch}$ competes with those from $f_\mathrm{DMI}$ to provide a nearly-collinear spin structure with a weak canting angle $\phi^\mathrm{WFM} = \cos^{-1}{\left(\mathbf{m}_1\cdot \mathbf{m}_2\right)}$. The introduction of an effective field due to $f_\mathrm{easy}$ puts the collinearity within the easy plane defined by the director $\hat{\mathbf{P}}$. This yields $\theta_1 = \theta_2 = 90^\circ$ where $\theta_\eta = \cos^{-1}{\left(\mathbf{m}_\eta \cdot \hat{\mathbf{P}}\right)}$. The weak effective field generated from $f_\mathrm{anis}$ splits the easy-plane degeneracy into a six-fold symmetry upon energy minimization. The `Materials` block,

!listing tutorial/BFO_P0A0_mRD.i
         block=Materials
         link=False
         language=python

seeds the values of these quantities for $f_\mathrm{mag}$ (note that they are set at the top of the input file for ease of use). This is because the MOOSE `Parser` reads in `${quantity}` from the out-of-block definition `quantity = something`. The definition for `alpha_long` corresponds to that of $\tilde\alpha_\parallel$ and we set it with a constant `ParsedFunction`. It could be defined the same way as the other materials coefficients but we leave this as an option for the user in case they want to give $\tilde\alpha_\parallel$ some functional dependence (on time or space for example). We use the value of $\gamma$ or $\gamma / \mu_0 M_s$ is set here where the units are given in nanometers microseconds over picograms.

The `Kernels` block is used to initialize the partial differential equation and computes the residual and jacobian contributions for each component of the sublattices,

!listing tutorial/BFO_P0A0_mRD.i
         block=Kernels
         link=False
         language=python

One can see the derivations of the residual and jacobian entries for the `Kernels` `Objects` by visiting the `Syntax` page. Some more options for this problem are located in the `Executioner` block (i.e. time step size `dt`),

!listing tutorial/BFO_P0A0_mRD.i
         block=Executioner
         link=False
         language=python

where we have used the [NewmarkBeta](https://mooseframework.inl.gov/source/timeintegrators/NewmarkBeta.html) time integration method. We find that this is most numerically stable for AFM ringdown problems. The simulation runs fairly quickly (400 seconds) on 4 processors using the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) distribution of MOOSE, stepping through a few thousands of time steps to ringdown $\mathbf{m}_\eta$. A visualization of the components of $\mathbf{m}_\eta$ from the `ParaView` filter `PlotGlobalVariablesOverTime` is provided in the below figure.

!media media/tut_ringdown.png style=display:block;margin:auto;width:60%; caption=Angular quantities $\phi^\mathrm{WFM}, \theta_1,$ and $\theta_2$ during ringdown. id=tut_AFM_ground

By setting the initial condition of $\mathbf{m}_\eta$ in different directions, one can obtain the below table for the eight-fold orientations of $\mathbf{P}\uparrow \mathbf{A}$.

!equation
\begin{array}{c|cccccccc}
\hline
\hline
 & & & & & & & & \\
 (\mathbf{P}\uparrow\uparrow\mathbf{A}): & [111] & [\bar{1}11] & [1\bar{1}\bar{1}] & [\bar{1}\bar{1}\bar{1}] & [1\bar{1}1] & [11\bar{1}] & [\bar{1}\bar{1}1] & [\bar{1}1\bar{1}] \\
 \hline
 & & & & & & & & \\
\mathbf{L} \simeq & [\bar{1}10] & [101] & [101] & [\bar{1}10] & [110] & [\bar{1}10] & [\bar{1}10] & [110] \\
& [\bar{1}01] & [110] & [110] & [\bar{1}01] & [011] & [011] & [011] & [011] \\
& [0\bar{1}1] & [01\bar{1}] & [01\bar{1}] & [0\bar{1}1] & [\bar{1}01] & [101] & [101] & [\bar{1}01] \\
& [1\bar{1}0] & [\bar{1}0\bar{1}] & [\bar{1}0\bar{1}] & [1\bar{1}0]  & [\bar{1}\bar{1}0] & [1\bar{1}0] & [1\bar{1}0] & [0\bar{1}\bar{1}] \\
& [10\bar{1}] & [\bar{1}\bar{1}0] & [\bar{1}\bar{1}0] & [10\bar{1}] & [0\bar{1}\bar{1}] & [0\bar{1}\bar{1}] & [0\bar{1}\bar{1}] & [0\bar{1}\bar{1}] \\
& [01\bar{1}] & [0\bar{1}1] & [0\bar{1}1] & [01\bar{1}] & [10\bar{1}] & [\bar{1}0\bar{1}] & [\bar{1}0\bar{1}] & [10\bar{1}] \\
\hline
& & & & & & & & & & & & & & & & \\
\mathbf{m} \simeq & [\bar{1}\bar{1}2] & [12\bar{1}] & [\bar{1}\bar{2}1] & [11\bar{2}] & [\bar{1}12] & [112] & [\bar{1}\bar{1}\bar{2}] & [1\bar{1}\bar{2}] \\
 & [1\bar{2}1] & [\bar{1}1\bar{2}] & [1\bar{1}2] & [\bar{1}2\bar{1}] & [\bar{2}\bar{1}1] & [2\bar{1}1] & [\bar{2}1\bar{1}] & [21\bar{1}] \\
 & [2\bar{1}\bar{1}] & [\bar{2}\bar{1}\bar{1}] & [211] & [\bar{2}11] & [\bar{1}\bar{2}\bar{1}] & [1\bar{2}\bar{1}] & [\bar{1}21] & [121]\\
 & [11\bar{2}] & [\bar{1}\bar{2}1] & [12\bar{1}] & [\bar{1}\bar{1}2] & [1\bar{1}\bar{2}] & [\bar{1}\bar{1}\bar{2}] & [112] & [\bar{1}12] \\
 & [\bar{1}2\bar{1}] & [1\bar{1}2] & [\bar{1}1\bar{2}] & [1\bar{2}1] & [21\bar{1}] & [\bar{2}1\bar{1}] & [2\bar{1}1] & [\bar{2}\bar{1}1] \\
 & [\bar{2}11] & [211] & [\bar{2}\bar{1}\bar{1}] & [2\bar{1}\bar{1}] & [121] & [\bar{1}21] & [1\bar{2}\bar{1}] & [\bar{1}\bar{2}\bar{1}] \\
\end{array}

Importantly, the nonzero value of $\mathbf{m}$ corresponds to a canted angle $\phi^\mathrm{WFM}$ of about $1.21^\circ$ which agrees well with the literature on BFO.


-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_results2.md next=MSCA_EU_Horizon2020_Results/horizon2020_ex2.md
