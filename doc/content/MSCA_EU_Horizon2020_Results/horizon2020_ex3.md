!config navigation breadcrumbs=False

# H2020-MSCA-IF-2019, Example: magnetic DWs

This page details how to obtain the magnetic textures in the presence of a FE domain boundary. We relax the Landau-Lifshitz-Bloch equation in the presence of large damping ($\alpha = 0.8$),

\begin{equation}\label{eqn:LLG_LLB}
  \begin{aligned}
    \frac{d \mathbf{m}_\eta}{dt} = -\frac{\gamma}{1+\alpha^2} \left(\mathbf{m}_\eta \times \mathbf{H}_\eta\right) - \frac{\gamma\alpha}{1+\alpha^2}\mathbf{m}_\eta\times\left(\mathbf{m}_\eta\times\mathbf{H}_\eta\right) + \frac{\gamma\tilde\alpha_\parallel}{(1+\alpha^2)} m_{\eta}^2 \left[ m_{\eta}^2 - 1 \right]\mathbf{m}_\eta,
  \end{aligned}
\end{equation}

where $\mathbf{H}_\eta = \mu_0^{-1} M_s^{-1} \delta f_\mathrm{mag} / \delta \mathbf{m}_\eta$ is the effective field. `The Kernels` block computes the relevant contributions from $\mathbf{H}_\eta$,

here

The parameters $\gamma$, $M_s$, and $\mu_0$ are the electron gyromagnetic, magnetization saturation density, and vacuum permeability. Note that $f_\mathrm{mag} = f_\mathrm{sp} + f_\mathrm{MP}$ is the total magnetic energy split between a component ($f_\mathrm{sp}$) that does not depend on the information of the structural order and a component that does ($f_\mathrm{MP}$).

Therefore, one needs to load in the $\mathbf{P}$ and $\mathbf{A}$ as an initial condition. For this, we use the previous example calculation but refined slightly (use `BFO_dwP1A1_100_ref.e` from the tutorials folder). We use refinement so that the magnetic solution to the LLG-LLB equation is numerically stable. We load with the Mesh block,

Next, we define the variables to solve for in the problem ($\mathbf{m}_\eta$). Note that $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$ are postprocessed as `AuxVariables`. We select an initial condition corresponding to a possible rotation of $\mathbf{L}$. However, this is not stable and the system evolves to minimize gradients in $\mathbf{L}$. Different


Relaxation of $\mathbf{m}_\eta$ with this mesh takes around 5260 seconds on 6 processors. We supply the reader with the output of this mesh file if needed.

The DW profile of $\mathbf{m}$ can be visualized in `ParaView` by using the `PlotOverLine` filter and could like the below figure,

!media media/tut_mag_DW.png style=display:block;margin:auto;width:60%; caption=Components of $\mathbf{m}$ across the 1/1 (100) DW. The net magnetization switches by $71^\circ$. id=tut_mag_DW

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_ex2.md next=MSCA_EU_Horizon2020_Results/horizon2020_ex4.md
