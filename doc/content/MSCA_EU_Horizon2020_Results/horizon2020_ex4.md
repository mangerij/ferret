!config navigation breadcrumbs=False

# H2020-MSCA-IF-2019, Example: ME switching

This page details how to obtain ME switching trajectories in FERRET. We use an electric field to switch $\mathbf{P}$ explicitly in a fully-time dependent simulation. The magnetization follows $\mathbf{P}$ since any change in $\mathbf{P}$ (and $\mathbf{A}$) corresponds to a drastic change in the magnetic energy, $f_\mathrm{MP}$, that couples to the structural distortions.

Our full approach involves solving the couple dynamic equation system,

\begin{equation}\label{eqn:TDLG}
  \begin{aligned}
    \frac{\partial \mathbf{P} }{\partial t} &= -\Gamma_P \frac{\delta f}{\delta \mathbf{P}}, \\
    \frac{\partial \mathbf{A} }{\partial t} &= -\Gamma_A \frac{\delta f}{\delta \mathbf{A}},
  \end{aligned}
\end{equation}

for the structural order ($\{\mathbf{P},\mathbf{A}\}$) along with,

\begin{equation}
  \begin{aligned}
    \frac{\partial \sigma_{ij}}{\partial x_j} = \frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + Q_{klmn} P_m P_n + R_{klmn} A_m A_n \right)\right] = 0,
  \end{aligned}
\end{equation}

at every time step, where $Q_{klmn}$ and $R_{klmn}$ are the electro- and rotostrictive coefficients that couple $\mathbf{P}$ and $\mathbf{A}$ to the elastic strain respectively. For the spins, we solve,

\begin{equation}\label{eqn:LLG}
  \begin{aligned}
    \frac{\partial \mathbf{M}_\eta}{\partial t} = -\left(\frac{\gamma}{1+\alpha^2}\right)\mathbf{M}_\eta\times \mathbf{H}_\eta - \frac{1}{M_s} \left(\frac{\gamma \alpha}{1+\alpha^2}\right) \mathbf{M}_\eta \times \left(\mathbf{M}_\eta \times \mathbf{H}_\eta\right).
  \end{aligned}
\end{equation}

Here, $\Gamma_P, \Gamma_A$ are a relaxation coefficients related to the time scales involved in the structural phase transition. The parameter $\gamma$ is the electron gyromagnetic ratio and $\mathbf{H}_\eta$ is the effective field acting on sublattice $\eta$. The coeffiicent $\alpha$ is a phenomenological damping constant which if made nonzero (and positive) drives the magnetic system to the ground state. We select a time-dependent electric field $\mathbf{E}$ of the form,

\begin{equation}\label{eqn:E}
  \begin{aligned}
    \mathbf{E}(\omega) = \langle 0,0, E_0 \rangle \sin{\left(\omega t\right)}.
  \end{aligned}
\end{equation}

which facilitates a reversal of the $P_z$ component. This means we expect the transition of $\mathbf{P}$ to go from $[111]$ to $[11\bar{1}]$ orientation. For this problem, we use the input file `BFO_P111_TO_P111b_switch_m1_a1.i` located in the tutorials subdirectory. The `Exodus` input that we use corresponds to that of the second simulation in Example 1 (a fully relaxed polar-magnetic solution). We load this output as an input via the `Mesh` block,

!listing tutorial/BFO_P111_TO_P111b_switch_m1_a1.i
         block=Mesh
         link=False
         language=python

with corresponding flags in the `Variables` block,

!listing tutorial/BFO_P111_TO_P111b_switch_m1_a1.i
         block=Variables
         link=False
         language=python

The `Kernels` are due to structural evolution (TDLGD) along with the ones from micromagnetic evolution (LLG-LLB). They are listed in the following lengthy block,

!listing tutorial/BFO_P111_TO_P111b_switch_m1_a1.i
         block=Kernels
         link=False
         language=python

For example, [ElectrostrictiveCouplingDispDerivative](source/kernels/ElectrostrictiveCouplingDispDerivative.md) and [RotostrictiveCouplingDispDerivative](source/kernels/RotostrictiveCouplingDispDerivative.md) correspond to the RHS of the mechanical equilibrium condition,

\begin{equation}
  \begin{aligned}
    \frac{\partial \sigma_{ij}}{\partial x_j} = \frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + Q_{klmn} P_m P_n + R_{klmn} A_m A_n \right)\right] = 0.
  \end{aligned}
\end{equation}

The `Kernel` [AFMEasyPlaneAnisotropySC](source/kernels/AFMEasyPlaneAnisotropySC.md) corresponds to the RHS of,

\begin{equation}\label{eqn:LLG2}
  \begin{aligned}
    \frac{\partial \mathbf{M}_\eta}{\partial t} = -\left(\frac{\gamma}{1+\alpha^2}\right)\mathbf{M}_\eta\times \mathbf{H}_\eta^\mathrm{easy-plane} - \frac{1}{M_s} \left(\frac{\gamma \alpha}{1+\alpha^2}\right) \mathbf{M}_\eta \times \left(\mathbf{M}_\eta \times \mathbf{H}_\eta^\mathrm{easy-plane}\right).
  \end{aligned}
\end{equation}

where $\mathbf{H}_\eta^\mathrm{easy-plane}$ is the effective field due to the free energy density term responsible for easy-plane magnetic anisotropy,

\begin{equation}\label{eqn:easy}
  \begin{aligned}
    f_\mathrm{easy-plane} = K_1 \left(\textbf{m}_\eta \cdot \hat{\mathbf{P}}\right)^2
  \end{aligned}
\end{equation}

which we compute as

\begin{equation}\label{eqn:easyeff}
  \begin{aligned}
    \mathbf{H}_\eta^\mathrm{easy-plane} = -\mu_0^{-1} M_s^{-1} \frac{\delta \left[ K_1 \left(\textbf{m}_\eta \cdot \hat{\mathbf{P}}\right)^2\right]}{\delta \mathbf{m}_\eta}
  \end{aligned}
\end{equation}

The suffix `SC` in the name of this object corresponds to a strongly-coupled situation. This means that for variational derivatives with respect to $\mathbf{P}$, this object also contributes non-zero residual contributions (i.e. in the time dependence of $\mathbf{P}$. We refer the reader to our `Syntax` page for a extensive list of the `Kernels` in FERRET regarding this problem. We also evaluate a number of postprocessed quantities in the `AuxKernels` block,

!listing tutorial/BFO_P111_TO_P111b_switch_m1_a1.i
         block=AuxKernels
         link=False
         language=python

For example, we compute $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$ with `VectorDiffOrSum`. We also compute other values such as, i.e., the angular quantities $\phi^\mathrm{WFM} = \cos^{-1}{\left(\mathbf{m}_1 \cdot \mathbf{m}_2\right)}$ and $\theta_\eta = \cos^{-1}{\left(\mathbf{m}_\eta \cdot \hat{\mathbf{P}}\right)}$. Note that we use the convention that $|\mathbf{L}| + |\mathbf{m}| = 2$ instead of unity. . The `Materials` block, assigns values to our coefficients,

!listing tutorial/BFO_P111_TO_P111b_switch_m1_a1.i
         block=Materials
         link=False
         language=python

where we have used units of `nanometers`, `microseconds`, `attocoulombs`, and `picograms`. This sets the time and length scales in this problem. The `Executioner` block, chooses flags for the time integration and numerical solve,

!listing tutorial/BFO_P111_TO_P111b_switch_m1_a1.i
         block=Executioner
         link=False
         language=python

A possible visualization of the output using ParaView is provided below,

!media media/tut_ME_sw.png style=display:block;margin:auto;width:67%; caption=Homogeneous switching at $\alpha = 0.003$ with Top: Neel vector switching and Bottom: net magnetization switching. These dynamics are acquired using a time-dependent electric field at frequency $\omega = 600$ MHz.  id=fig_ME_sw

The wall clock time for this problem is 62.37 seconds on 6 processors using the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) distribution of MOOSE. Note that other switching trajectories can be obtained by switching out the initial `Exodus` file - for example, choosing different six-fold ${\mathbf{L},\mathbf{m}}$ orientation or by selecting a different $\mathbf{E}$ orientation.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_ex3.md next=MSCA_EU_Horizon2020_Results/horizon2020_ex5.md
