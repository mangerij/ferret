# Micromagnetic theory

!alert construction title=Documentation in-progress
This section requires some work before it will be finalized online. Please contact the developers if you need assistance with this aspect of the module.

The micromagnetic approach to modeling magnetic systems is similar to the phase field method in the sense that a continuum description of an order parameter is provided. In ferromagnets, below the critical temperature, the spins order forming a nonzero and net spin magnetization $\mathbf{M}$. Additionally, the spins see a specific symmetry characterized by the presence (or absence) of magnetocrystalline anisotropy. This drives the tendency to form ordered domains but competitions from the ferromagnetic exchange or the demagnetizing fields (internal or external) can favor other spin textures such as vortices, spin cycloids, or skymrions. To predict the magnetic ground states (and/or their dynamics), the Landau-Lifshitz-Gilbert (LLG) equation is used,

\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{M}}{\partial t} = -\gamma \mathbf{M}\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M} \times \mathbf{M}\times \mathbf{H}
  \end{aligned}
\end{equation}

where $\gamma$ is the electron gyromagnetic ratio and $\alpha$ is a phenomenological damping coefficient. The vector $\mathbf{H}$ is the effective field which arises from different physics of the magnet (exchange, anisotropy, spin-torque transfer...). The second term in the governing equation drives the system towards its energy minima (or system equilibrium).

In contrast to ferroelectrics where the magnitude of the electric dipole moment $\mathbf{P}$ can stretch or shrink, the magnetization $\mathbf{m} = \mathbf{M} / M_s$ stays on the unit sphere during the LLG time evolution.

There are two possible numerical resources in FERRET/MOOSE to use to enforce this constraint. One possibility is to renormalize $\mathbf{m}$ at the conclusion of every time step in the evolution of Eq. (1). Another option is to evolve the Landau-Liftshitz-Bloch (LLB) equation

\begin{equation}\label{LLG-LLB}
  \begin{aligned}
    \frac{\partial \mathbf{M}}{\partial t} = -\gamma \mathbf{M}\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M} \times \mathbf{M}\times \mathbf{H} - \frac{\gamma \alpha_\mathrm{LLB}}{1+\alpha^2} \left(m^2 - 1\right) \mathbf{m}
  \end{aligned}
\end{equation}

in which the latter term provides a longitudinal restoring force to $\mathbf{m}$ mapping it back onto the unit sphere. We have demonstrated that both methods are essentially equivalent (they also can be used together).

If the material is an antiferromagnet, then the spin system forms two sublattices $(\mathbf{M}_1, \mathbf{M}_2)$ whose interaction is mediated by a short-range (typically nearest-neighbor) exchange. To simulate the micromagnetics of an antiferromagnet, we have the two sublattice LLG equation,

\begin{equation}
  \begin{aligned}
    \frac{\partial \mathbf{M}_\eta}{\partial t} = -\gamma \mathbf{M}_\eta\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M}_\eta \times \mathbf{M}_\eta \times \mathbf{H},
  \end{aligned}
\end{equation}

where $\eta = 1, 2$. As mentioned before, one can choose to use the LLB equation or let the modulus be renormalized at every time step. In a collinear antiferromagnet, the total net magnetization, $\mathbf{M} = \mathbf{M}_1 + \mathbf{M}_2$, is zero indicating that the spin sublattices perfectly compensating each other. 
