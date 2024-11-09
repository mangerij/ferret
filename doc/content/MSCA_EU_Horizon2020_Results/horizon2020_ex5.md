!config navigation breadcrumbs=False

# H2020-MSCA-IF-2019, Example: spin wave transport

This page details how to obtain spin wave transport across the multiferroic domain boundary. In order to perform this calculation, one needs the information from two other examples (first relax the FE DW, then get the relaxed magnetic DW in the presence of the FE domain boundary). These are handled in Example 2 and 3 of the MSCA tutorials respectively.

The initial condition is a fixed $\{\mathbf{P},\mathbf{A}\}$ configuration with a corresponding relaxed $\{\mathbf{L},\mathbf{m}\}$ ground state. We load in the `Exodus` file from the previous example in the `Mesh` block,

!listing tutorial/BFOspw_P111A111-P111bA111b_m1_H111_alp0a0f7DW.i
         block=Mesh
         link=False
         language=python

and the `Variables` block,

!listing tutorial/BFOspw_P111A111-P111bA111b_m1_H111_alp0a0f7DW.i
         block=Variables
         link=False
         language=python

which uses the `initial_from_file_var` with the `LATEST` for `initial_from_file_timestep` option. The values of $\mathbf{L}$ and $\mathbf{m}$ are calculated in the `AuxKernels`. We then use the `ParsedFunction` system to set up an applied perturbation to generate spin waves in the $\pm\mathbf{k}$ directions,

!listing tutorial/BFOspw_P111A111-P111bA111b_m1_H111_alp0a0f7DW.i
         block=Functions
         link=False
         language=python

with `pulse1`, `pulse2`, and `pulse3` being the $x$, $y$, and $z$ components of $\mathbf{H}_\mathrm{appl}$ defined as,

\begin{equation}\label{eqn:perturb}
  \begin{aligned}
  \mathbf{H}_\mathrm{appl} &= H_0 \, \mathrm{sinc}[k_0 (x-x_0)] \, \mathrm{sinc}[\omega_0 (t-t_0)] \, e^{-p_0(x-x_0)^2} \,\hat{\mathbf{h}}
  \end{aligned}
\end{equation}

where field amplitude $H_0$, excitation location $x_0$, gaussian intensity profile parameter $p_0 = 0.16$ $\mathrm{nm}^{-2}$, and $k_0 = 10$ $\mathrm{nm}^{-1}$ control the perturbation distribution in spacetime. The director $\hat{\mathbf{h}}$ orients the magnetic field with respect to $\mathbf{m}$. Finally, we cut-off the pulse at $t_0 = 1$ ps and excite the spin waves at a frequency $\omega_0 = 500$ GHz. The `Kernels` for this problem are set as,

!listing tutorial/BFOspw_P111A111-P111bA111b_m1_H111_alp0a0f7DW.i
         block=Kernels
         link=False
         language=python

which is the same as the previous magnetic DW simulation of Example 3 with the exception of `AFMInteractionCartLLHConst` which adds the residual and jacobian contributions due to the Zeeman interaction,

\begin{equation}\label{eqn:zeeman}
  \begin{aligned}
   f_\mathrm{Zeeman} = -\mathbf{m}\cdot\mathbf{H}_\mathrm{appl}
  \end{aligned}
\end{equation}

for a constant $\mathbf{H}_\mathrm{appl}$ field (no dependence on $\Phi_\mathrm{H}$). The `Executioner` options are set,

!listing tutorial/BFOspw_P111A111-P111bA111b_m1_H111_alp0a0f7DW.i
         block=Executioner
         link=False
         language=python

with the `NewmarkBeta` time integration method and an upper limit of $dt = 50$ fs. A visualization of the spin wave traveling through the DW region is provided below from `ParaView`,

!media media/sw_ex_page.mp4 style=display:block;margin:auto;width:67%; caption=Animation of spin wave transport through the 1/1 (100) FE domain boundary. The plotted variables are the components of the weak net magnetization $m_k$. id=fig_sw_ex_page

This problem has a wall-clock time of 276.46 seconds on 6 processors using the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) distribution of MOOSE. Different spin waves can be generate by changing the definitions of the `ParsedFunctions`.



-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez-González](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_ex4.md
