# Background

!alert construction title=Documentation in-progress
This section will have information that will be available when the review process is complete for the main publication encompassing this work. Please contact the developers if you need assistance with this aspect of the module.

Within the scope of the [Project SCALES - 897614](https://cordis.europa.eu/project/id/897614), the idea was to develop a new continuum approach for simulation of multiferroic compounds.  Multiferroics are a class of materials which display both an electric and magnetic ordering below their critical phase transition. Specifically, we use the perovskite $\mathrm{BiFeO}_3$ (BFO) as our material workhorse. BFO displays both ferroelectricity and antiferromagnetism at room temperature, leading to a host of possible applications from beyond-CMOS logic gates, tunneling magnetoresistant spintronic valves, THz radiation emitters, enhanced piezoelectric elements, ultrafast acoustic modulators, to linear electrooptical components.

To study these materials at the device-relevant scale in this project, we proposed to couple the phase field method for ferroelectrics with the micromagnetic approach. The idea was to have, self-consistently, the same time and length scale for both electric and magnetic order unlocking the magnetoelectric coupling in both static or dynamic configurations (i.e. switching or spin-wave transport) in a scalable open-source code.

In order to disseminate results to the public in the framework of the H2020-MSCA-IF-2019 call, this additional page was added to the FERRET website. We will discuss background and the model, summarize our key findings, and provide examples of how to reproduce the results obtained during the funding period and described in our article [!cite](Mangeri2023).

Stay tuned for an update!

# Properties of $\mathrm{BiFeO}_3$ and model

We consider a zero temperature limit free energy density functional defined as a sum of Landau-Devonshire energy density from the structural distortions of the lattice ($f_\mathrm{lat}$), the magnetic energy density due to the spin subsystem ($f_\mathrm{sp}$), and the ME coupling between the order parameters ($f_\mathrm{ME}$) in single crystal BFO.

\begin{equation}\label{eqn:en_sum}
  \begin{aligned}
    f &= f_\mathrm{lat}(\mathbf{P},\mathbf{A},\bm{\varepsilon}) + f_\mathrm{sp}(\mathbf{L},\mathbf{m}) + f_\mathrm{ME}(\mathbf{L},\mathbf{m}, \mathbf{P},\mathbf{A})
  \end{aligned}
\end{equation}

In the continuum description, some formal definitions of the order parameters are needed. The electric polarization $\mathbf{P}$ is connected to the displacement of the $\mathrm{Bi}^{3+}$ and $\mathrm{Fe}^{3+}$ relative to the oxygen anions along the pseudocubic $\[111\]$ or equivalent directions. The vector $\mathbf{A}$ describes the rotations of the $\mathrm{FeO}_6$ cages where the antiphase behavior in adjacent unit cells is implicitly assumed. See Fig. \ref{OPsPA}.

!media media/OPsPA.png style=display:block;margin:auto;width:60%; caption=Structural distortion order parameters $\mathbf{P}$ and $\mathbf{A}$ corresponding to the electric dipole moment and antiphase octahedral rotations respectively.  id=OPsPA

The spontaneous homogeneous strain that arises below the phase transition is the rank two tensor $\bm{\varepsilon}$ with symmetric components $\varepsilon_{ij} = \varepsilon_{ji}$,

\begin{equation}\label{eqn:strain}
  \begin{aligned}
    \varepsilon_{ij} &= \frac{1}{2} \left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)
  \end{aligned}
\end{equation}

The variable $u_i$ is the component of the elastic displacement vector $\mathbf{u}$ accompanies the ferroelectric phase transition in this material (sometimes called spontaneous or ferroelastic strain).

For the spin subsystem, BFO is an antiferromagnet with anti-aligned spins at adjacent Fe sites (G-type) leading to two distinct sublattices $\mathbf{m}_1$ and $\mathbf{m}_2$. The quantity $\mathbf{L}$ is the AFM N\'{e}el vector which we define as $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$. The total magnetic moment $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$ which is the weak nonvanishing magnetization that arises due to an antisymmetric Dzhaloshinksii-Moriya interaction (DMI) - see [!cite](Ederer2005) and related work.

!media media/OPsM.png style=display:block;margin:auto;width:60%; caption=Atomistic depiction of the magnetic sublattices $\mathbf{m}_\eta$ for $\eta = 1,2$ demonstrating canting due to DMI in the easy-plane.  id=OPsM

The quantities $\mathbf{L}$ and $\mathbf{m}$ are defined such that $|\mathbf{L}| + |\mathbf{m}| = 1$ with, in general, $|\mathbf{L}| \gg |\mathbf{m}|$ reflecting the presence of a strong AFM coupling between the sublattices but with a weak noncollinearity in $\mathbf{m}_1$ and $\mathbf{m}_2$. The total weak magnetization is $\mathbf{M} = M_s \mathbf{m}$ where $M_s$ is the saturation magnetization density of the Fe sublattice ($\approx 4.0 \mu$ B/Fe) [!cite](Dixit2015).

A detailed description of our model is shared in 

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded in 2021-2023 at the Luxembourg Institute of Science and Technology under the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) call.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;
