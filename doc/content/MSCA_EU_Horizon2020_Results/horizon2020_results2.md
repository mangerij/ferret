!config navigation breadcrumbs=False

# Results (Model)

The goals of the project were to first develop a model for $\mathrm{BiFeO}_3$ BFO that couples the electric polarization, oxygen octahedral cage antiphase tilt structure, and the magnetization suitable for evaluating the multiferroic domain topology and dynamics. The first task that we had was to parameterize the energy penalty for the formation of domain walls (DWs) against predictions from density functional theory (DFT). Once this was accomplished, we were able to capture the inhomogeneous magnetic order in the presence of the structural DWs. This concludes our investigation of the ground states and development of the model. Next we used it in two different simulations corresponding to fully-coupled dynamical switching along with spin wave transport through the multiferroic domain boundary. We summarize each of these pillars of the project below. We refer the reader to the previous page for some definitions of the properties of BFO.

# Homogeneous structural orders

This section details our results of homogeneous solutions of BFO for the structural and spin order. We first concern ourselves with the structural order that results at the conclusion of the evolution of the time dependent Landau-Ginzburg equations,

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

where $Q_{ijkl}$ and $R_{ijkl}$ are the electrostrictive and rotostrictive coefficient tensors with $\mathbf{P}$ and $\mathbf{A}$ the order parameters associated with the spontaneous polarization and antiphase tilting of the oxygen octahedral cages. The mechanical equilibrium equation is assumed to be satisfied for every time step which is a reasonable approximation for the elastic strains that arise during domain evolution. The total free energy density, $f$ corresponds to the fourth-order coupled potential parameterized by the DFT work detailed [!cite](Fedorova2022). We set $\Gamma_P = \Gamma_A$ to unity for these calculations since we are only interested in the final state. Initial conditions are set such that $\mathbf{P}$ and $\mathbf{A}$ are small quantities near zero and the boundary conditions (periodic) are both chosen to produce a homogeneous solution with $\mathbf{P}$ and $\mathbf{A}$ parallel ($\mathbf{P}\uparrow\uparrow\mathbf{A}$). The magnitudes of the structural order parameters are $P_s = |\mathbf{P}| = 0.945 \,\,\mathrm{C}/\mathrm{m}^2$  and $A_s = |\mathbf{A}| = 13.398^\circ$. The spontaneous $homogeneous$ normal and shear strain values are listed as rows with $\varepsilon_n = 1.308\times10^{-2}$ and $\varepsilon_s = 2.95 \times 10^{-3}$ respectively. The eight-fold domain variant symmetry is found in the pseudocubic reference frame corresponding to the below table.

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

To reproduce these results, we refer the reader [here](MSCA_EU_Horizon2020_Results/horizon2020_ex1.md) for more details on the input file and the documentation to reproduce these results.

# Antiferromagnetic ringdown (homogeneous spin state)

To find the possible magnetic states in this material, we evolve the normalized two-sublattice Landau-Lifshitz-Bloch (LLB) equation,

\begin{equation}\label{eqn:LLG_LLB}
  \begin{aligned}
    \frac{d \mathbf{m}_\eta}{dt} = -\frac{\gamma}{1+\alpha^2} \left(\mathbf{m}_\eta \times \mathbf{H}_\eta\right) - \frac{\gamma\alpha}{1+\alpha^2}\mathbf{m}_\eta\times\left(\mathbf{m}_\eta\times\mathbf{H}_\eta\right) + \frac{\gamma\tilde\alpha_\parallel}{(1+\alpha^2)} m_{\eta}^2 \left[ m_{\eta}^2 - 1 \right]\mathbf{m}_\eta,
  \end{aligned}
\end{equation}

where $\mathbf{P}$ and $\mathbf{A}$ are held fixed corresponding to one of the possible eight-fold degenerate domain directions. This approximation for evaluation of ground states is suitable since the magnetic system does not appreciably influence the relaxation dynamics of the polarization since the work required to move ions (i.e. those relative displacements corresponding to $\mathbf{P}$ and $\mathbf{A}$) is much larger than those to reorient the spins on the Fe sites. The constant $\gamma$ corresponds to the electron gyromagnetic factor equal to $2.2101\times10^5$ rad. m A${}^{-1}$ s ${}^{-1}$, $\alpha$ the Gilbert damping, and $\tilde\alpha_\parallel$ the longitudinal damping constant from the LLB approximation. The effective fields are defined as $\mathbf{H}_\eta = - \mu_0^{-1} M_s^{-1} \delta f / \delta \mathbf{m}_\eta$ with $\mu_0$ the permeability of vacuum. The saturation magnetization density of the BFO sublattices ($\eta = 1,2$) is $M_s = 4.0$ $\mu$B/Fe [!cite](Dixit2015). We set $\alpha = 0.01$ and $\tilde\alpha_\parallel = 10^2$ and ringdown the magnetic system to an energy minimum.

The below figure details the numerical ringdown behavior of the sublattices $\mathbf{m}_\eta$. Here, we also look for homogeneous spins (no domain walls or texture in the $\mathbf{m}_\eta$ field).

!media media/ringdown_abc.png style=display:block;margin:auto;width:100%; caption=Ringdown time dependence of Left: components of selected $\mathbf{m}_\eta$. Center: easy-plane angle $\theta_\eta = \cos^{-1}{(\mathbf{m}_\eta\cdot\hat{\mathbf{P}})}$ and Right: canted moment angle $\phi^\mathrm{WFM} = \cos^{-1}{(\mathbf{m}_1\cdot \mathbf{m}_2)}$.  id=fig_ringdown_abc

Importantly, this analysis allowed us to parameterize the DMI coupling at our level of theory, which leads to the non-vanishing moment of $\mathbf{m} = (\mathbf{m}_1 + \mathbf{m}_2)/2$. By tracking the angle of the canted moment $\phi^\mathrm{WFM} = \cos^{-1}{\left(\mathbf{m}_1 \cdot \mathbf{m}_2 \right)}$ as a function of the strength of the DMI coupling, $D_0$, we can identify good agreement ($\phi^\mathrm{WFM} \approx 1.22^\circ$) with the available literature corresponding to $\mathbf{M} = M_s\mathbf{m} = 0.03 \mu$B/f.u..

By selecting initial conditions corresponding to the six-fold possible orientations of the spin sublattice system, the below table can be obtained yielding 48 different orientations of the the Néel order $\mathbf{L} = \left(\mathbf{m}_1 - \mathbf{m}_2\right) / 2$ and total magnetization $\mathbf{m}$ dependent on the 8 possible orientations of $(\mathbf{P}\uparrow\uparrow\mathbf{A})$.

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

To reproduce these results, we refer the reader to our examples documentation [here](MSCA_EU_Horizon2020_Results/horizon2020_ex1.md) for details on the simulation (input/output).

# Structural domain walls

The lowest-order invariants of gradients in the structural fields are,

\begin{equation}\label{eqn:gradP}
  \begin{aligned}
    f_\mathrm{\nabla P} =\frac{G_{11}}{2}  \left( P_{x,x}^2 + P_{y,y}^2 + P_{z,z}^2 \right) +  G_{12}  \left(P_{x,x} P_{y,y} + P_{y,y} P_{z,z} + P_{x,x} P_{z,z} \right) + \frac{G_{44}}{2} \left[\left(P_{x,y} + P_{y,x} \right)^2+ \left(P_{y,z} + P_{z,y} \right)^2 + \left(P_{x,z} + P_{z,x}\right)^2\right]
  \end{aligned}
\end{equation}

and

\begin{equation}\label{eqn:gradA}
  \begin{aligned}
    f_\mathrm{\nabla A} = \frac{H_{11}}{2}   \left( A_{x,x}^2 + A_{y,y}^2 + A_{z,z}^2 \right) +  H_{12}  \left(A_{x,x} A_{y,y} + A_{y,y} A_{z,z} + A_{x,x} A_{z,z} \right) + \frac{H_{44}}{2} \left[\left(A_{x,y} + A_{y,x} \right)^2+ \left(A_{y,z} + A_{z,y} \right)^2 + \left(A_{x,z} + A_{z,x}\right)^2\right].
  \end{aligned}
\end{equation}

A comma in the subscript denotes a partial derivative with respect to a specified spatial direction. In order to study the domain wall topology involving spatial variations of $\mathbf{P}$, $\mathbf{A}$, and strain, a good parameter set estimate of the gradient coefficients $(G_{11}, H_{11}, ...)$ is needed. To achieve this, we consult DFT calculations reported in [!cite](Dieguez2013). It was shown that an assortment of different metastable (along with low energy) domain walls form an energetic hierarchy. This was advantageous for our work, as we were allowed to seperate the computation (fitting) of various $(G_{11}, H_{11}, ...)$ due to only certain terms in the gradient energy density expressions being identified as primary contributions to the DW energy.

For example, consider the 2/1 (100) DW which is a commonly observed domain boundary observed in experiments. In this notation, it is indicated that, for the $2/1$ DW, two components of $\mathbf{P}$ and one component of $\mathbf{A}$ vary across the boundary whose plane normal is (100) whereis for the $3/0$ DW, $\mathbf{P}$ undergoes a full reversal where $\mathbf{A}$ is approximately unchanged across the (110)-oriented boundary plane. We label the pairs of the domains characterizing the DW as $\mathbf{P}^\mathrm{I}/\mathbf{A}^\mathrm{I}$ and $\mathbf{P}^\mathrm{II}/\mathbf{A}^\mathrm{II}$ in the Table below.

After seeding a $sin(x)$ profile for the initial state of $\mathbf{P}$ and $\mathbf{A}$ corresponding to the expected DW configuration, a typical profile of the order parameters and strain can be obtained as shown below for the $2/1 [100]$ DW when relaxing the time-dependent Landau-Ginzburg equations - Eqs. (\ref{eqn:TDLGD}).

!media media/tut_21_DW.png style=display:block;margin:auto;width:47%; caption=Different components of $\mathbf{P}$ and $\mathbf{A}$ across the $2/1 [100]$ domain wall.  id=fig_dw_21_prof

We compute the DW energy with $F_\mathrm{DW} = \left(F - F_0\right)/\left(N S\right)$ where $F_0$ is the corresponding monodomain energy density from the fourth-order potential integrated over the computational volume. The energy $F$ is computed from the solution that contains the DW profile with the number of DWs in the simulation box being $N$ and $S$ the surface area of the DW plane. Therefore, $F_\mathrm{DW}$ can be computed as a function of $G_{ij}, H_{ij}$ for each wall. This allowed us to separate different contributions of the gradient coefficients. We extend this type of analysis iteratively throughout the possible DWs so that we can converge our set of coefficients yielding reasonable $F_\mathrm{DW}$ values comparable to DFT; importantly, capturing the energy hierarchy predicted (by [!cite](Dieguez2013)) for metastable walls.

The energy hierarchy as compared to results from DFT is shown in the below Table

!equation
\begin{array}{ccccc}
\hline
\hline
 & & & & & \\
\textbf{P}^\mathrm{I}/\mathbf{A}^\mathrm{I} & \mathrm{Type} & \mathrm{DW} & \textbf{P}^\mathrm{II}/\mathbf{A}^\mathrm{II} & F_\mathrm{DW}^{(\mathrm{DFT})} & F_\mathrm{DW}^{(\mathrm{FEM})}  \\
 & & & & & \\
 \hline
[111]/[111] & 0/0 & - & [111]/[111] & - & -  \\
[111]/[111] & 0/3 & (100) & [111]/[\bar{1}\bar{1}\bar{1}] & 227 & 293  \\
[111]/[111] & 1/1 & (100) & [11\bar{1}]/[11\bar{1}] & 151 & 162 \\
[111]/[111] & 1/2 & (100) & [11\bar{1}]/[\bar{1}\bar{1}1] & 147 & 159  \\
[111]/[111] & 2/1 & (100) & [1\bar{1}\bar{1}]/[\bar{1}11] & 62 & 60  \\
[111]/[111] & 2/2 & (100) & [1\bar{1}\bar{1}]/[1\bar{1}\bar{1}] & 319 & 314 \\
[1\bar{1}1]/[1\bar{1}1] & 3/0 & (110) & [\bar{1}1\bar{1}]/[1\bar{1}1] & 74 & 78 \\
[1\bar{1}1]/[1\bar{1}1] & 3/3 & (110) & [\bar{1}1\bar{1}]/[\bar{1}1\bar{1}] & 255 & 263 \\
\end{array}

To reproduce the profiles of a 2/1 (100) DW, we refer the reader to our examples documentation [here](MSCA_EU_Horizon2020_Results/horizon2020_ex2.md). Other domain walls are easily obtained by switching out the `ICs`.

# Domain walls in magnetization

After we obtained the structural DW profile, the two sublattice LLG-LLB equation,

\begin{equation}\label{eqn:LLG}
  \begin{aligned}
    \frac{\partial \mathbf{M}_\eta}{\partial t} = -\frac{\gamma}{1+\alpha^2} \mathbf{M}_\eta\times \mathbf{H}_\eta - \frac{1}{M_s} \left(\frac{\gamma \alpha}{1+\alpha^2}\right) \mathbf{M}_\eta \times \mathbf{M}_\eta \times \mathbf{H}_\eta,
  \end{aligned}
\end{equation}

is relaxed with large $\alpha$ to find the ground states of $\mathbf{M}_\eta$ in the presence of the FE domain boundary. In this equation, $\gamma$ is the electron gyromagnetic ratio and,

\begin{equation}\label{eqn:eff}
  \begin{aligned}
    \mathbf{H}_\eta = - \frac{1}{\mu_0 M_s}\frac{\delta f}{\delta \mathbf{m}_\eta},
  \end{aligned}
\end{equation

is the effective field which is calculated locally by variational derivatives of the free energy density. We include a non-local exchange free energy density contribution,

\begin{equation}
  \begin{aligned}
    f_\mathrm{\nabla L} = A_e \left[\left(\nabla L_x\right)^2 + \left(\nabla L_y\right)^2 + \left(\nabla L_z\right)^2\right].
  \end{aligned}
\end{equation}

where $A_e = 3 \times 10^{-7}$ ergs/cm as proposed in the work of [!cite](Agbelele2017). The resulting magnetization texture is shown below.

!media media/tut_11_DW.png style=display:block;margin:auto;width:47%; caption=Different components of $\mathbf{m}$ across the $1/1 [100]$ domain wall.  id=fig_dw_11_mag

Here, we see in black that the component that does not switch has a slight bending indicating a small rotation. Since the $\mathbf{A}$ vector switches by $71^\circ$, this forces the weak magnetization moment $\mathbf{m} = (\mathbf{m}_1 + \mathbf{m}_2) / 2$ to also rotate by $71^\circ$.

To generate these results, we refer the reader to our examples documentation [here](MSCA_EU_Horizon2020_Results/horizon2020_ex3.md).  Other domain walls are easily obtained by switching out the `ICs` for $\mathbf{m}_\eta$ or choosing different starting DWs of the $\mathbf{P}$ and $\mathbf{A}$ subsystem.

# Magnetoelectric switching simulations

In AFM spintronics, it is needed to find ways to control the spin order using external stimuli. However, when working with noncollinear AFM BFO, the magnetization is weak, making it difficult to use a uniform applied magnetic field through the Zeeman interaction. To overcome this challenge, an electric field can be used to manipulate and control the magnetic texture since BFO has an intrinsic electric dipole moment. The potential benefits of electric field control of magnetism have been studied for some time, and while low-frequency deterministic switching of $\mathbf{P}$ has been demonstrated in the work of [!cite](Heron2014), the dynamic processes of the coupled polar-spin order at higher frequencies still require further research. This study focuses on the use of modeling to investigate ME switching, which involves using an electric field to switch $\left\{\mathbf{L},\mathbf{m}\right\}$.

We consider a fully-dynamical simulation where all system variables $\mathbf{P},\mathbf{A},\mathbf{u},\Phi_\mathrm{E},\mathbf{m}_1$, and $\mathbf{m}_2$ depend on time. The dynamics of the AFM order are in general very fast ($100$s of GHz to THz regime), therefore, we introduce a numerical constraint on the time steps for the evolution of Eq. (\ref{eqn:LLG_LLB}) of dt $ < 0.1$ ps. To switch the $z$ component of $\mathbf{P}$, we choose our electric field to be $\mathbf{E}(\omega) = \langle 0,0, E_0 \sin{\left(\omega t\right)}\rangle$. We select an $\mathbf{E}$ frequency of $\omega = 600$ MHz. The field is abruptly turned off after $\mathbf{P}$ has switched in order to facilitate only one switching event for analysis. The initial state is homogeneous $\mathbf{P}\uparrow\uparrow\mathbf{A}$ along $[111]$ with $\mathbf{L}||[\bar{1}01]$ and $\mathbf{m}||[1\bar{1}1]$ as one of the possibilities listed in the above Table of possible magnetic ground states.

In our paper, we investigated the dependence of the switching trajectories on the Gilbert damping $\alpha$ which for magnetic insulators (and BFO) is postulated to be on the order of $10^{-3}$. Shown in the below figure, the choice of $\alpha$ strongly influences the final states of $\mathbf{L}$ despite having very similar ringdown profiles in the vicinity of the $\mathbf{P}$ switching event.

!media media/tut_SW.png style=display:block;margin:auto;width:67%; caption=Homogeneous switching with Left: $\alpha = 0.003$ and Right: $\alpha = 0.01$. The dynamics of the Neel vector $\mathbf{L}$ are acquired using a time-dependent electric field at frequency $\omega = 600$ MHz.  id=fig_sw

Not shown here is the time dependence of the weak magnetization $\mathbf{m}$ which also switches in both cases (see article).

To reproduce these results, we refer the reader to our examples documentation [here](MSCA_EU_Horizon2020_Results/horizon2020_ex4.md). Other switching trajectories are readily obtained by changing the `ICs` for $\mathbf{m}_\eta$ and orientations of $\mathbf{E}$.

# Spin wave transport across the multiferroic domain boundary

Spintronics involves creating, managing, and detecting spin packets, which are important for information processing (see [!cite](Hirohata2020)). In AFMs, spin precessional processes occur at low energies and very high frequencies, offering advantages over standard CMOS technology(see [!cite](Jungwirth2016)). Understanding wave transmission and reflection is important for optimizing spin wave transport in these systems. A recent study from [!cite](Parsonet2022) has shown that thermal magnon transport in BFO can be controlled using electric fields. The study found that the $109^\circ$ FE DWs act as a barrier to spin transport over a length-scale of several hundred nanometers, hindering the detected magnon signal that is useful for the device.

We investigated this particular scenario mesoscopically by leveraging our model described in the previous sections. We consider the two of the commonly observed DWs in BFO experiments, the $109^\circ$ 2/1 and $71^\circ$ 1/1 (100)-oriented walls. There is a large relative difference in energies between the lattice and spin contribution  (i.e. $|f_\mathrm{latt}| \gg |f_\mathrm{sp} + f_\mathrm{MP}|$). This suggests that any application of an external magnetic field $\mathbf{H}_\mathrm{appl}$ should not appreciably influence the $\mathbf{P}$ and $\mathbf{A}$ subsystem. Therefore, we fix (in time) these order parameters in this section. We couple $\mathbf{H}_\mathrm{appl}$ to act on the weak moment through the Zeeman free energy density,

\begin{equation}\label{eqn:Zeeman}
  \begin{aligned}
    f_\mathrm{Zeeman} = -\mathbf{m}\cdot\mathbf{H}_\mathrm{appl}
  \end{aligned}
\end{equation}

and add it to the total free energy of the spin configuration. To generate spin waves, we consider gaussian perturbations generated by the following field form (see [!cite](Gruszecki2015)),

\begin{equation}\label{eqn:perturb}
  \begin{aligned}
  \mathbf{H}_\mathrm{appl} &= H_0 \, \mathrm{sinc}[k_0 (x-x_0)] \, e^{-p_0(x-x_0)^2} \mathrm{sinc}[\omega_0 (t-t_0)] \,\mathbf{h}
  \end{aligned}
\end{equation}

where field amplitude $H_0 = 184$ kOe, excitation location $x_0$, gaussian intensity profile parameter $p_0 = 0.16$ $\mathrm{nm}^{-2}$, and $k_0 = 10$ $\mathrm{nm}^{-1}$ control the perturbation distribution in spacetime. The director $\mathbf{h}$ orients the magnetic field with respect to $\mathbf{m}$. Finally, we cut-off the pulse at $t_0 = 1$ ps and excite the spin waves at a frequency $\omega_0$. Eq. (\ref{eqn:LLG_LLB}) is evolved with $\alpha = 0$ and Eq. (\ref{eqn:perturb}).

We enforce periodicity in our computational volume along the $x, y, z$ for the $\mathbf{m}_1$ and $\mathbf{m}_2$ variables. The time-integration of Eq.~(\ref{eqn:LLG_LLB}) is set for $dt < 2$ fs time steps to ensure numerical convergence for the fast AFM dynamics in the system. We verify that our calculations are in the linear limit by adjusting the $H_0$ and determining that the perturbed amplitudes of $\mathbf{m}_\eta$ scale linearly. Finally, we monitor the system total free energy $F_\mathrm{sp} + F_\mathrm{ME}$ and $|\mathbf{m}_\eta|$ (via the LLB term) and verify that they are both constant for all time in our $\alpha = 0$ simulation.

In the below video, we plot the excess energy as as a function of the arclength (perpendicular) to the DW in nanometers. The animation shows the time dependence of the spin wave package as it hits the DW located at 22 nm.

!media media/sw_ex.mp4 style=display:block;margin:auto;width:67%; caption=Spin wave transport through the 2/1 (100) FE domain boundary. The color map on the Gylphs is set for $m_y = m_z$ which changes sign across the DW. The right image is a plot along the arclength $[100]$ of the components of $\mathbf{m}$. id=fig_sw_ex

To reproduce these results, we refer the reader to our examples documentation [here](MSCA_EU_Horizon2020_Results/horizon2020_ex5.md). Other spin wave transport calculations can be set up by changing the initial perturbing field or changing the `ICs` for $\{\mathbf{P},\mathbf{A},\mathbf{m}_\eta\}$ corresponding to any of the DWs described above.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_results1.md next=MSCA_EU_Horizon2020_Results/horizon2020_ex1.md
