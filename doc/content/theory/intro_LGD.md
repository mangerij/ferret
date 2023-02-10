# Phase field method of ferroelectrics

!alert construction title=Documentation in-progress
This section requires some work before it will be finalized online. Please contact the developers if you need assistance with this aspect of the module

!media media/fig1_doublewell.png style=display:block;margin:auto;width:50%; caption=Left; conventional double well energy surface of a ferroelectric below $T_C$. Right; Nonzero $P_s$ as a function of temperature.  id=fig1_doublewell

The ferroelectric (FE) phase field method utilizes a thermodynamic approach to describe the FE phase transition. Within this theoretical description, these materials are characterized by the onset of a primary order parameter that becomes nonzero below a critical Curie temperature $T_C$. In the case of canonical FEs such as $\mathrm{PbTiO}_3$ or $\mathrm{BaTiO}_3$, the primary order parameter is the electric dipole moment. In both of these compounds, below $T_C$, the material symmetry is lowered from cubic to tetragonal. The ionic displacements (the relative off-centering the Ti atom) give rise, spontaneously, to a net nonzero polarization or $\mathbf{P}_s$.

To describe this transition, the quasi-dynamical Landau-Gilbert-Devonshire (LGD) equation can be evolved to search for ground states of the system,

\begin{equation}\label{LGD}
  \begin{aligned}
    \frac{\partial \mathbf{P}}{\partial t} = - \Gamma_P \frac{\delta F}{\delta \mathbf{P}}.
  \end{aligned}
\end{equation}

Here, the system energy is given by $F$ whose minimum is when $\mathbf{P} = \mathbf{P}_s$. The action of evolving Eq. (1) is to follow the energy down towards the minimum by starting from a sufficiently close estimate of the cubic state with $\mathbf{P}(\mathbf{r}) \approx 0$. Since the FE material has more than one energy minima corresponding to the tetragonal (six-fold) symmetry (see above double well for $F$, domains will form characterized by equal energy configurations of $\pm\mathbf{P}$. The strength of the FE phase field method is that it allows the identification of arbitrary orientations of $\mathbf{P}$ in the continuum limit (agnostic of length scales typically restrictive of atomistic simulations). The domain patterns can be predicted given a relatively simple set of inputs. Within the FERRET/MOOSE ecosystem, the finite element method is employed to discretize Eq. (1) onto regular or irregular grids. The use of the latter allows for one to search for domain patterns in arbitrary geometries (i.e. nanoparticles or curved surfaces).

In thin film samples, these domain orientations depend strongly on temperature or other external applied fields (i.e. electric or mechanical in origin) which can also be included in the calculations. The coupling to the electric fields (internal and/or external) are provided by the Poisson equation,

\begin{equation}\label{Poisson}
  \begin{aligned}
    \nabla^2 \Phi_\mathrm{E} = \nabla \cdot \mathbf{P}.
  \end{aligned}
\end{equation}

where the electric field is defined in the usual way $\mathbf{E} = - \nabla \Phi_\mathrm{E}$.

If the system has a strong dependence on elastic fields (as is the case for most ferroelectrics), then mechanical equilibrium can be sought be solving

\begin{equation}\label{stressdiv}
  \begin{aligned}
    \frac{\sigma_{ij}}{\partial x_j} = 0,
  \end{aligned}
\end{equation}

where $\sigma_{ij}$ is the total stress tensor and $\partial / \partial x_j$ denotes spatial derivatives in the $x_j$ direction. Eqs (2) and (3) can be solved at every time step in the evolution of Eq. (1). In the phase field approach, the decomposition of $F$ contains different contributions due to different physics,

\begin{equation}\label{totalenergy}
  \begin{aligned}
    f = f_\mathrm{bulk} + f_\mathrm{grad} + f_\mathrm{elastic} + f_\mathrm{elec} + f_\mathrm{electrostrict} + ...,
  \end{aligned}
\end{equation}

with $f_\mathrm{bulk}$ the bulk double well potential density, $f_\mathrm{grad}$ the gradient energy density penalty for formation of domain walls, $f_\mathrm{elastic}$ linear elastic energy density, $f_\mathrm{elec}$ interaction energy density for the inclusion of an internal or external field, $f_\mathrm{electrostric}$ the electrostrictive energy density coupling between the dipole moment and the strain fields. All of these terms can be expanded up to arbitrary order depending on the material symmetry and need for accuracy of predicted spontaneous order parameter values.

For the gradient energy density, one typically uses the lowest-order Lifshitz invariants as described in [!cite](Cao1990) and [!cite](Hlinka2006),

\begin{equation}\label{lifshitz}
  \begin{aligned}
    f_\mathrm{grad} &=\frac{1}{2} G_{11}  \left\{ \left(\frac{\partial P_x}{\partial x} \right)^2 + \left(\frac{\partial P_y}{\partial y} \right)^2+\left(\frac{\partial P_z}{\partial z} \right)^2 \right\} +  G_{12}  \left\{\left(\frac{\partial P_x}{\partial x} \right)\left(\frac{\partial P_y}{\partial y} \right) + \left(\frac{\partial P_y}{\partial y} \right)\left(\frac{\partial P_z}{\partial z} \right) + \left(\frac{\partial P_x}{\partial x} \right)\left(\frac{\partial P_z}{\partial z} \right)\right\} \\
  \end{aligned}
\end{equation}

which requires knowledge of $G_{ij}$. Higher order terms are possible although not generally neccessary to describe the domain wall structure in bulk or thin film.


For canonical perovskites, there is a good understanding of the values of these parameters.

!media media/time_evol.png style=display:block;margin:auto;width:50%; caption=Time evolution of the ferroelectric energy for a small block of material - competition arises from the formation of interfacial regions (domain walls) and the bulk energy which prefers aligned polar states.  id=time_evol
