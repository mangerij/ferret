# H2020-MSCA-IF-2019

!alert construction title=Documentation in-progress
This section will have information that will be available when the review process is complete for the main publication encompassing this work. Please contact the developers if you need assistance with this aspect of the module.

In the framework of the Marie-Curie individual fellowship - [H2020-MSCA-IF-2019](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) - call, this additional page was added to the FERRET website in order to disseminate results of [Project SCALES - 897614](https://cordis.europa.eu/project/id/897614) to the general public and the technical audience. We will discuss background and the modeling effort, summarize our key findings, and provide examples of how to reproduce representative results obtained during the funding period which are described in detail within our article [!cite](Mangeri2023) currently on the arXiv preprint server and under review.

# Background

The project proposal submitted in 2019 encapsulated the main idea of developing a new continuum approach for simulations of multiferroic compounds. Multiferroics are a class of materials which display both an electric and magnetic ordering below their critical phase transition temperature. Specifically, we use the perovskite $\mathrm{BiFeO}_3$ (BFO) as our material workhorse. BFO displays both ferroelectricity and antiferromagnetism at room temperature, leading to a host of possible applications from beyond-CMOS logic gates, tunneling magnetoresistant spintronic valves, THz radiation emitters, enhanced piezoelectric elements, ultrafast acoustic modulators, to linear electrooptical components. We refer the reader to INSERT for an excellent review of the properties of BFO along with its technological applications.

To study these materials at the device-relevant scale, we proposed to couple the phase field method for ferroelectrics with the micromagnetic approach. Both of these methodologies utilize a "coarse-grained" description of order parameters. Both approaches are phenomenological in origin and have been shown to describe well the theoretical picture of ferroelectric and ferromagnetic materials respectively.  Typically, this is challenging because the energy of the magnetic order tends to be many orders of magnitude lower than the structural distortions. Incidentally, this leads to different time and length scales which accompany the structural or magnetic phase transitions. The idea was to have, self-consistently on the same time and length scale, a model for both electric and magnetic order unlocking the magnetoelectric properties in both static or dynamic configurations. This allowed us to investigate some applications of the model corresponding to magnetoelectric switching and also spin-wave transport across the multiferroic domain boundaries.

We leveraged the Multiphysics Object Oriented Simulation Environment (MOOSE) framework which is open-source software convenient for rapid model development with advanced features. [MOOSE](https://mooseframework.inl.gov/) is developed and maintained at Idaho National Laboratory in the United States. We implement our models of BFO within FERRET (this website), an open-source add-on module (part of the MOOSE ecosystem of applications). See the [landing page](index.md) for more information on FERRET.

Stay tuned for an update!

# Properties of $\mathrm{BiFeO}_3$ and coupled model

We consider a zero temperature limit free energy density functional defined as a sum of Landau-type energy density from the structural distortions of the lattice ($f_\mathrm{latt}$), the magnetic energy density due to the nominally collinear spin subsystem ($f_\mathrm{sp}$), and a coupling between the electric and magnetic order parameters ($f_\mathrm{MP}$) in single crystal BFO.

\begin{equation}\label{eqn:en_sum}
  \begin{aligned}
    f &= f_\mathrm{latt}(\mathbf{P},\mathbf{A},\varepsilon) + f_\mathrm{sp}(\mathbf{L},\mathbf{m}) + f_\mathrm{MP}(\mathbf{L},\mathbf{m}, \mathbf{P},\mathbf{A})
  \end{aligned}
\end{equation}

In the continuum description, some formal definitions of the order parameters are needed. The electric polarization $\mathbf{P}$ is connected to the displacement of the $\mathrm{Bi}^{3+}$ and $\mathrm{Fe}^{3+}$ cations relative to the oxygen anions along the pseudocubic $[111]$ or equivalent directions. The vector $\mathbf{A}$ describes the rotations of the $\mathrm{FeO}_6$ cages about the polar axis $\mathbf{P}$ where the antiphase behavior in first-neighboring unit cells is implicitly assumed. See Fig. 1 which shows both of these structural distortions.

!media media/OPsPA.png style=display:block;margin:auto;width:60%; caption=Structural distortion order parameters $\mathbf{P}$ and $\mathbf{A}$ corresponding to the electric dipole moment and antiphase oxygen octahedral cage rotations about the polar axis respectively.  id=OPsPA

A homogeneous strain arises below the structural phase transition temperature which is a rank two tensor $\varepsilon$ with symmetric components $\varepsilon_{ij} = \varepsilon_{ji}$,

\begin{equation}\label{eqn:strain}
  \begin{aligned}
    \varepsilon_{ij} &= \frac{1}{2} \left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)
  \end{aligned}
\end{equation}

The variable $u_i$ is the component of the elastic displacement vector $\mathbf{u}$ which accompanies the ferroelectric phase transition in this material (sometimes called spontaneous or ferroelastic strain).

For the spin subsystem, BFO is an antiferromagnet with anti-aligned spins at first-neighboring Fe sites (G-type) leading to two distinct sublattices $\mathbf{m}_1$ and $\mathbf{m}_2$. The quantity $\mathbf{L}$ is the AFM N\'{e}el vector which we define as $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$. The magnetic orientation is of the easy-plane variety, whose plane normal is the polar director. We propose that the presence of the antisymmetric Dzhaloshinksii-Moriya interaction (DMI) introduces a free energy density term of the form,

\begin{equation}\label{eqn:dmi}
  \begin{aligned}
    f_\mathrm{DMI} &= D_0 \mathbf{A} \cdot \mathbf{m} \times \mathbf{L}
  \end{aligned}
\end{equation}

coupled to the direction of the antiphase tilts $\mathbf{A}$ which acts to cant the sublattices slightly in the plane. We refer the reader to [!cite](Ederer2005) and related work for more information on this phenomena in BFO. Therefore, due to the DMI, a weak nonvanishing magnetic moment is present in the system and can be computed as $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$.

!media media/OPsM.png style=display:block;margin:auto;width:60%; caption=Atomistic depiction of the magnetic sublattices $\mathbf{m}_\eta$ for $\eta = 1,2$ demonstrating canting due to DMI in the easy-plane (cyan circles).  id=OPsM

A schematic of the spin order is shown in Fig. 2 (left) easy-plane anisotropy of $\mathbf{m}_\eta$ relative to $\mathbf{P}$ direction and (right) DMI-induced canting in the easy-plane (shown as cyan circles). The quantities $\mathbf{L}$ and $\mathbf{m}$ are defined such that $|\mathbf{L}| + |\mathbf{m}| = 1$ with, in general, $|\mathbf{L}| \gg |\mathbf{m}|$ reflecting the presence of a strong AFM coupling between the sublattices but with a weak noncollinearity in $\mathbf{m}_1$ and $\mathbf{m}_2$. The total weak magnetization is $\mathbf{M} = M_s \mathbf{m}$ where $M_s$ is the saturation magnetization density of the Fe sublattice ($\approx 4.0 \mu$ B/Fe) [!cite](Dixit2015).

We solve the couple dynamic equation system,

\begin{equation}\label{eqn:TDLG}
  \begin{aligned}
    \frac{\partial \mathbf{P} }{\partial t} &= -\Gamma_P \frac{\delta f}{\delta \mathbf{P}}, \\
    \frac{\partial \mathbf{A} }{\partial t} &= -\Gamma_A \frac{\delta f}{\delta \mathbf{A}},
  \end{aligned}
\end{equation}

for the structural order ($\{\mathbf{P},\mathbf{A}\}$) and

\begin{equation}\label{eqn:LLG}
  \begin{aligned}
    \frac{\partial \mathbf{M}_\eta}{\partial t} = -\gamma \mathbf{M}_\eta\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M}_\eta \times \mathbf{M}_\eta \times \mathbf{H},
  \end{aligned}
\end{equation}

for the spin system. Here, $\Gamma_P, \Gamma_A$ are a relaxation coefficients related to the time scales involved in the phase transition. The coeffiicent $\alpha$ is a phenomenological damping constant which if made nonzero (and positive) drives the magnetic system to the ground state.

A detailed description of our model is shared in the preprint on arXiv at [!cite](Mangeri2023).

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;
